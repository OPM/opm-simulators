#ifndef DUNE_RANDOMPERMEABILITY_HH
#define DUNE_RANDOMPERMEABILITY_HH

#include<iostream>
#include<vector>
#include<set>
#include<map>
#include<stdio.h>
#include<stdlib.h>

#include <dune/disc/functions/p0function.hh>

namespace Dune
{
	/*! \brief providing the absolute permeability.
	 *
	 *  The class Permeability is derived from the template argument BV which usually 
	 *  respresents a block vector. The values for the permeability should be already set by 
	 *  the constructor. 
	 */
	template<class G>
	class RandomPermeability {
	    template<int dim>
	    struct ElementLayout
	    {
	      bool contains (Dune::GeometryType gt)
	      {
		return gt.dim() == dim;
	      }
	    }; 
	  
		enum{n = G::dimension};
		typedef typename G::ctype DT;
		typedef LeafP0Function<G,DT,1> PermType;
		typedef BlockVector<FieldVector<DT,1> > RepresentationType;
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::LeafGridView GV;
	    typedef typename GV::IndexSet IS;
	    typedef Dune::MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> EM;
	
	public:
	/*! \brief Constructor.
	 *
	 *  \param size number of degrees of freedom   
	 *  \param grid pointer to grid 
	 *  \param mapper pointer to mapper 
	 *
	 *  The constructor already sets the entries to 
	 *  permeability values which are described by the functor \a permFunc, 
	 *  as for example given by PermeabilityBall or RandomPermeability. If the bool 
	 *  \a permFunc.random is set to true, \a permFunc is expected to set all entries in one call. 
	 *  Otherwise, a traversal over the cells is done, and \a permFunc should return the 
	 *  permeability at the cell center. 
	 */
		RandomPermeability(const G& g, const char* name = "permeab.dat", const bool create = true)
		: grid(g), perm(g), permloc(0), 
		  createNew(create), fileName(name), elementmapper(g, g.leafIndexSet())
		{
			typedef typename GV::template Codim<0>::Iterator Iterator;
			    	
			const GV& gridview(grid.leafView());
		    Iterator eendit = gridview.template end<0>();
			
		    char* pwd(getenv("PWD"));
		    char startCommand[220];
		    strcpy(startCommand, "find ");
		    strcat(startCommand, pwd);
		    strcat(startCommand, "/");
		    char systemCommand[220];

		    char simsetdir[221];
		    int foundSimSet = 0; 
		    int k = 0;
		    while (!foundSimSet && k++ < 5) {
		      strcpy(systemCommand, startCommand);
		      strcat(systemCommand, " -type f -name simset -exec dirname {} \\; >simsetloc.txt");
		      system(systemCommand);
		      std::ifstream simsetloc("simsetloc.txt");
		      simsetloc.seekg (0, std::ios::end);
		      int length = simsetloc.tellg();
		      simsetloc.seekg (0, std::ios::beg);
		      if (length > 0) {
			foundSimSet = 1;
			simsetloc.getline(simsetdir, 220);
		      }
		      simsetloc.close();
		      strcat(startCommand, "../");
		    }
		    
		    if (createNew)
		    {
			    // SIMSET creates random permeabilities for given coordinates, so the coordinates of the center of gravity of each element 
			    // are written to a file 'SIMKOR'
			    // open output stream for simset output file name 
			    char namefileName[100];
			    strcpy(namefileName, simsetdir);
			    strcat(namefileName, "/SIMNAM");
			    std::ofstream namefile(namefileName);
			    // Choose simset output filename
			    namefile << fileName << std::endl;
			    namefile.close();
			    // open output stream for simset input file
			    char outfileName[100];
			    strcpy(outfileName, simsetdir);
			    strcat(outfileName, "/SIMKOR");
			    std::ofstream outfile(outfileName);
			    for (Iterator it = gridview.template begin<0>(); it != eendit; ++it)
			    {
			    	Dune::GeometryType gt = it->geometry().type(); 
		
					const Dune::FieldVector<DT,n>& 
					  local = Dune::ReferenceElements<DT,n>::general(gt).position(0,0);
					
					// get global coordinate of cell center
					Dune::FieldVector<DT,n> global = it->geometry().global(local);
					
				    outfile << global[0] << "\t" << global[1] << std::endl;
			    }
				outfile.close();
				strcpy(systemCommand, "cd ");
				strcat(systemCommand, simsetdir);
				strcat(systemCommand, "; ./simset; cd $OLDPWD");
				system(systemCommand);
		    }
		    
			// open input stream for simset output file
		        char concd[100];
			strcpy (concd, simsetdir);
			strcat(concd, "/");
			std::ifstream infile(strcat(concd, fileName));
			std::cout << "Read permeability data from " << concd << std::endl;
			for (Iterator it = gridview.template begin<0>(); it != eendit; ++it)
			{
				int indexi = elementmapper.map(*it);
				double dummy1, dummy2, permi;
				char zeile [221];
				infile.getline(zeile, 220);
				std::istringstream ist(zeile);
				ist >> dummy1 >> dummy2 >> permi;
				(*perm)[indexi] = pow(10.0, permi);
		    }
		    infile.close();
		}
		
		//! return const reference to permeability vector
		const RepresentationType& operator* () const
		{
		  return (*perm);
		}
	
		//! return reference to permeability vector
		RepresentationType& operator* ()
		{
		  return (*perm);
		}
	
		Dune::FieldMatrix<DT,n,n>& K (const Entity& e) 
		{
			int elemId = elementmapper.map(e);
			DT permE = (*perm)[elemId];

			for (int i = 0; i < n; i++)
				permloc[i][i] = permE;

			return permloc;
		}
		
	  void vtkout (const char* name, const G& grid) const 
	  {
	    Dune::VTKWriter<G, typename G::LeafGridView> 
	      vtkwriter(grid.leafView());
	    vtkwriter.addCellData(*perm, "absolute permeability");
	    int size = (*perm).size();
	    RepresentationType logPerm(size);
	    for (int i = 0; i < size; i++) 
	      logPerm[i] = log10((*perm)[i]);
	    vtkwriter.addCellData(logPerm, "logarithm of permeability");
	    vtkwriter.write(name, Dune::VTKOptions::ascii);		
	  }
	
	private:
		const G& grid; 
		PermType perm;
		Dune::FieldMatrix<DT,n,n> permloc;
		const bool createNew;
		const char* fileName;
		EM elementmapper;
	};
	
	
	
	
	
	/*! \brief providing the absolute permeability for cells on given level and their children.
	 *
	 *  Unlike the class RandomPermeability 
	 *  which provides the permeability for the leaf grid, in LevelRandomPermeability the
	 *  permeability field is provided on a given grid level \f$ l \f$. The permeability
	 *  of the level-\f$ l \f$ elements is also inherited to their children. 
	 */
	template<class G>
	class LevelRandomPermeability {
	    template<int dim>
	    struct ElementLayout
	    {
	      bool contains (Dune::GeometryType gt)
	      {
		return gt.dim() == dim;
	      }
	    }; 
	  
		enum{n = G::dimension};
		typedef typename G::ctype DT;
		typedef LevelP0Function<G,DT,1> PermType;
		typedef BlockVector<FieldVector<DT,1> > RepresentationType;
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::Traits::template Codim<0>::EntityPointer EntityPointer;
		typedef typename G::LevelGridView GV;
	    typedef typename GV::IndexSet IS;
	    typedef Dune::MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> EM;
	
	public:
	/*! \brief Constructor.
	 *
	 *  \param g a grid objectof type G  
	 *  \param lev the level on which the permeability is to be provided 
	 *  \param name the name of the file in the simset-directory in which the permeabilities are to be stored. 
	 *  \param create set true if new field shall be created, set false if permeabilities shall be read from specified file.
	 * 
	 *  The constructor already sets the entries to 
	 *  permeability values which are described by the functor \a permFunc, 
	 *  as for example given by PermeabilityBall or RandomPermeability. If the bool 
	 *  \a permFunc.random is set to true, \a permFunc is expected to set all entries in one call. 
	 *  Otherwise, a traversal over the cells is done, and \a permFunc should return the 
	 *  permeability at the cell center. 
	 */
		LevelRandomPermeability(const G& g, const int lev, const char* name = "permeab.dat", const bool create = true)
		: grid(g), perm(g,lev), permloc(0), level_(lev),
		  createNew(create), fileName(name), elementmapper(g, g.levelIndexSet(lev))
		{
			if (lev > g.maxLevel() ) DUNE_THROW(Dune::Exception,"Level specified for permeability data is higher than maximum grid level!");
			typedef typename GV::template Codim<0>::Iterator Iterator;
			    	
			const GV& gridview(grid.levelView(level()));
		    Iterator eendit = gridview.template end<0>();
			
		    char* pwd(getenv("PWD"));
		    char startCommand[220];
		    strcpy(startCommand, "find ");
		    strcat(startCommand, pwd);
		    strcat(startCommand, "/");
		    char systemCommand[220];

		    char simsetdir[221];
		    int foundSimSet = 0; 
		    int k = 0;
		    while (!foundSimSet && k++ < 5) {
		      strcpy(systemCommand, startCommand);
		      strcat(systemCommand, " -type f -name simset -exec dirname {} \\; >simsetloc.txt");
		      system(systemCommand);
		      std::ifstream simsetloc("simsetloc.txt");
		      simsetloc.seekg (0, std::ios::end);
		      int length = simsetloc.tellg();
		      simsetloc.seekg (0, std::ios::beg);
		      if (length > 0) {
			foundSimSet = 1;
			simsetloc.getline(simsetdir, 220);
		      }
		      simsetloc.close();
		      strcat(startCommand, "../");
		    }
		    
		    if (createNew)
		    {
			    // SIMSET creates random permeabilities for given coordinates, so the coordinates of the center of gravity of each element 
			    // are written to a file 'SIMKOR'
			    // open output stream for simset output file name 
			    char namefileName[100];
			    strcpy(namefileName, simsetdir);
			    strcat(namefileName, "/SIMNAM");
			    std::ofstream namefile(namefileName);
			    // Choose simset output filename
			    namefile << fileName << std::endl;
			    namefile.close();
			    // open output stream for simset input file
			    char outfileName[100];
			    strcpy(outfileName, simsetdir);
			    strcat(outfileName, "/SIMKOR");
			    std::ofstream outfile(outfileName);
			    for (Iterator it = gridview.template begin<0>(); it != eendit; ++it)
			    {
			    	Dune::GeometryType gt = it->geometry().type(); 
		
					const Dune::FieldVector<DT,n>& 
					  local = Dune::ReferenceElements<DT,n>::general(gt).position(0,0);
					
					// get global coordinate of cell center
					Dune::FieldVector<DT,n> global = it->geometry().global(local);
					
				    outfile << global[0] << "\t" << global[1] << std::endl;
			    }
				outfile.close();
				strcpy(systemCommand, "cd ");
				strcat(systemCommand, simsetdir);
				strcat(systemCommand, "; ./simset; cd $OLDPWD");
				system(systemCommand);
		    }
		    
			// open input stream for simset output file
		        char concd[100];
			strcpy (concd, simsetdir);
			strcat(concd, "/");
			std::ifstream infile(strcat(concd, fileName));
			std::cout << "Read permeability data from " << concd << std::endl;
			for (Iterator it = gridview.template begin<0>(); it != eendit; ++it)
			{
				int indexi = elementmapper.map(*it);
				double dummy1, dummy2, permi;
				char zeile [221];
				infile.getline(zeile, 220);
				std::istringstream ist(zeile);
				ist >> dummy1 >> dummy2 >> permi;
				(*perm)[indexi] = pow(10.0, permi);
		    }
		    infile.close();
		}
		
		//! return const reference to permeability vector
		const RepresentationType& operator* () const
		{
		  return (*perm);
		}
	
		//! return reference to permeability vector
		RepresentationType& operator* ()
		{
		  return (*perm);
		}
	
		//! \brief return reference to permeability tensor of specified cell.
		/** \param e cell of level\f$ l \f$ or higher
		 * 
		 */
		Dune::FieldMatrix<DT,n,n>& K (const Entity& e) 
		{
			int le = e.level();
			int elemId;
			if (le < level_) DUNE_THROW(Dune::Exception, "Level of element lower than level of permeability discretisation, permeability not uniquely defined"); 
			else if (le > level_)
			{
				EntityPointer f = e.father();
				le = f->level();
				while (le > level_)
				{
					f = f->father();
					le = f->level();
				}
				elemId = elementmapper.map(*f);
			}
			else elemId = elementmapper.map(e);
			DT permE = (*perm)[elemId];

			for (int i = 0; i < n; i++)
				permloc[i][i] = permE;

			return permloc;
		}
		
	  void vtkout (const char* name, const G& grid) const 
	  {
	    Dune::VTKWriter<G, typename G::LevelGridView> 
	      vtkwriter(grid.levelView(level_));
	    int size = (*perm).size();
	    vtkwriter.addCellData(*perm, "absolute permeability");
	    RepresentationType logPerm(size);
	    for (int i = 0; i < size; i++) 
	      logPerm[i] = log10((*perm)[i]);
	    vtkwriter.addCellData(logPerm, "logarithm of permeability");
	    vtkwriter.write(name, Dune::VTKOptions::ascii);		
	  }
	
	  int level()
	  {
		  return level_;
	  }
	  
	private:
		const G& grid;
		const int level_;
		PermType perm;
		Dune::FieldMatrix<DT,n,n> permloc;
		const bool createNew;
		const char* fileName;
		EM elementmapper;
	};
}

#endif


