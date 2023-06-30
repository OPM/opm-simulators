#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <limits>
#include <algorithm>

#define BOOST_TEST_MODULE TestTuningXXXMBE
#include <boost/test/unit_test.hpp>

struct Column : public std::vector<std::string> {
		Column(const std::string& name_, const int size = 0, const int num_rows_estimate = 1000) : 
				std::vector<std::string>(size), name(name_)
		{
				this->reserve(num_rows_estimate);
		}

		// Return vector of double values, invalid elements set to NaN
		std::vector<double> dvalues() const {
				std::vector<double> vec;
				vec.reserve(this->size());
				
				const auto& conv_func = [](const std::string& strval) {
						double dval;
						try {
								dval = std::stod(strval);
						} catch (std::invalid_argument& exc) {
								dval = std::numeric_limits<double>::quiet_NaN();
						}
						return dval;
				};

				std::transform(this->cbegin(), this->cend(), std::back_inserter(vec), conv_func);
				return vec;
		}

		// Return vector of double values values, invalid elements set to std::numeric_limits<int>::min()
		std::vector<int> ivalues() const {
				std::vector<int> vec;
				vec.reserve(this->size());
				
				const auto& conv_func = [](const std::string& strval) {
						int ival;
						try {
								ival = std::stoi(strval);
						} catch (std::invalid_argument& exc) {
								ival = std::numeric_limits<int>::min();
						}
						return ival;
				};

				std::transform(this->cbegin(), this->cend(), std::back_inserter(vec), conv_func);
				return vec;
		}

		std::string name;
};


struct ColumnData {
		ColumnData(const std::string& file_name, const int num_columns_estimate=20) {
				raw_columns.reserve(num_columns_estimate);
				load_file(file_name);
		}

		void load_file(const std::string& file_name) {
				// Open file and read first line with column names
				std::ifstream ifs(file_name);
				std::string line, colname;

				std::getline(ifs, line);
				std::istringstream iss(line);
				while (iss >> colname) {
						column_names.push_back(colname);
						raw_columns.emplace_back(colname);
						columns[colname] = &(raw_columns.back());
				}
				const int num_columns = column_names.size();

				// Read remaining lines into std::string vectors
				int lineno = 1;
				while (std::getline(ifs, line)) {
						iss.str(line); iss.clear();
						int i=0;
						while (iss >> colname && i < num_columns) {
								raw_columns[i].push_back(colname);
								++i;
						}
						if (i >= num_columns && iss >> colname) {
								std::cout << "Warning: Ignoring extra column(s) on line " << lineno << std::endl;
						}
						++lineno;
				}
		}

		// Get data vectors of different types
		std::vector<double> get_dvector(const std::string& colname) const { return columns.at(colname)->dvalues(); }
		std::vector<int> get_ivector(const std::string& colname) const { return columns.at(colname)->ivalues(); }
		// Default is to return double values
		std::vector<double> operator[](const std::string& colname) const { return columns.at(colname)->dvalues(); }
		
		std::vector<std::string> column_names;
		std::vector<Column> raw_columns;
		std::map<std::string, Column*> columns;
};

BOOST_AUTO_TEST_CASE(CheckMassBalanceWithinXXXMBE)
{
		//std::string case_name(boost::unit_test::framework::master_test_suite().argv[1]);
		std::string case_name("01_TUNING_XXXMBE");
		std::string file_name = case_name + ".INFOITER";

		ColumnData data(file_name);
		auto rstep = data.get_ivector("ReportStep");
		auto tstep = data.get_ivector("TimeStep");
		auto mbo = data["MB_Oil"];
		auto mbw = data["MB_Water"];
		auto mbg = data["MB_Gas"];

		const int num_reports = 1 + *std::max_element(rstep.begin(), rstep.end());
		std::vector<double> max_mb;
		max_mb.reserve(num_reports);


		// Find the maximum mass balance error at each converged time step for each report step.. 
		const int nrows = rstep.size();
		int rcur = 0;
		int tcur = 0;
		double max_mb_step = std::numeric_limits<double>::min();
		for (int i=0; i<(nrows-1); ++i) {
				if (tcur != tstep[i+1] || rcur != rstep[i+1]) {
						max_mb_step = std::max({mbo[i], mbw[i], mbg[i], max_mb_step});
						tcur = tstep[i+1];
				}
				if (rcur != rstep[i+1]) {
						max_mb.push_back(max_mb_step);
						max_mb_step = std::numeric_limits<double>::min();
						rcur = rstep[i+1];
				}			
		}
		max_mb.push_back( std::max({mbo.back(), mbw.back(), mbg.back(), max_mb_step}) );

		BOOST_TEST_MESSAGE("---------------------------------------------------------------------------");
		BOOST_TEST_MESSAGE("Found the following converged max mass balance error (per report step):");
		for (auto& val : max_mb)
				BOOST_TEST_MESSAGE(val);
		BOOST_TEST_MESSAGE("---------------------------------------------------------------------------");


		BOOST_CHECK( max_mb[0] < 1.0e-6 );
		BOOST_CHECK( max_mb[1] < 1.0e-8 );
		BOOST_CHECK( max_mb[2] < 1.0e-10 );
  
}

