// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::VtkMultiWriter
 */
#ifndef EWOMS_VTK_MULTI_WRITER_HH
#define EWOMS_VTK_MULTI_WRITER_HH

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/bvector.hh>

#include <opm/material/common/Valgrind.hpp>

#include <opm/models/io/baseoutputwriter.hh>
#include <opm/models/io/vtkscalarfunction.hh>
#include <opm/models/io/vtktensorfunction.hh>
#include <opm/models/io/vtkvectorfunction.hh>

#include <opm/models/parallel/tasklets.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace Opm {

/*!
 * \brief Simplifies writing multi-file VTK datasets.
 *
 * This class automatically keeps the meta file up to date and
 * simplifies writing datasets consisting of multiple files. (i.e.
 * multiple time steps or grid refinements within a time step.)
 */
template <class GridView, Dune::VTK::OutputType vtkFormat>
class VtkMultiWriter : public BaseOutputWriter
{
    class WriteDataTasklet : public TaskletInterface
    {
    public:
        explicit WriteDataTasklet(VtkMultiWriter& multiWriter)
            : multiWriter_(multiWriter)
        {}

        void run() override final
        {
            std::string fileName;
            // write the actual data as vtu or vtp (plus the pieces file in the parallel case)
            if (multiWriter_.commSize_ > 1) {
                fileName = multiWriter_.curWriter_->pwrite(/*name=*/multiWriter_.curOutFileName_,
                                                           /*path=*/multiWriter_.outputDir_,
                                                           /*extendPath=*/"",
                                                           vtkFormat);
            }
            else {
                fileName = multiWriter_.curWriter_->write(/*name=*/multiWriter_.outputDir_ + "/" +
                                                          multiWriter_.curOutFileName_,
                                                          vtkFormat);
            }

            // determine name to write into the multi-file for the
            // current time step
            // The file names in the pvd file are relative, the path should therefore be stripped.
            const std::filesystem::path fullPath{fileName};
            const std::string localFileName = fullPath.filename();
            multiWriter_.multiFile_.precision(16);
            multiWriter_.multiFile_ << "   <DataSet timestep=\"" << multiWriter_.curTime_ << "\" file=\""
                                    << localFileName << "\"/>\n";
        }

    private:
        VtkMultiWriter& multiWriter_;
    };

    enum { dim = GridView::dimension };

    using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

public:
    using Scalar = BaseOutputWriter::Scalar;
    using Tensor = BaseOutputWriter::Tensor;
    using Vector = BaseOutputWriter::Vector;
    using ScalarBuffer = BaseOutputWriter::ScalarBuffer;
    using TensorBuffer = BaseOutputWriter::TensorBuffer;
    using VectorBuffer = BaseOutputWriter::VectorBuffer;

    using VtkWriter = Dune::VTKWriter<GridView>;

    VtkMultiWriter(bool asyncWriting,
                   const GridView& gridView,
                   const std::string& outputDir,
                   const std::string& simName = "",
                   const std::string& multiFileName = "")
        : gridView_(gridView)
        , elementMapper_(gridView, Dune::mcmgElementLayout())
        , vertexMapper_(gridView, Dune::mcmgVertexLayout())
        , outputDir_(outputDir.empty() ? "." : outputDir)
        , simName_(simName.empty() ? "sim" : simName)
        , multiFileName_(multiFileName.empty() ? outputDir_ + "/" + simName_ + ".pvd" : multiFileName)
        , commSize_(gridView.comm().size())
        , commRank_(gridView.comm().rank())
        , curWriterNum_(0)
        , taskletRunner_(/*numThreads=*/asyncWriting ? 1 : 0)
    {}

    ~VtkMultiWriter() override
    {
        taskletRunner_.barrier();
        releaseBuffers_();
        finishMultiFile_();

        if (commRank_ == 0) {
            multiFile_.close();
        }
    }

    /*!
     * \brief Returns the number of the current VTK file.
     */
    int curWriterNum() const
    { return curWriterNum_; }

    /*!
     * \brief Updates the internal data structures after mesh
     *        refinement.
     *
     * If the grid changes between two calls of beginWrite(), this
     * method _must_ be called before the second beginWrite()!
     */
    void gridChanged()
    {
        elementMapper_.update(gridView_);
        vertexMapper_.update(gridView_);
    }

    /*!
     * \brief Called whenever a new time step must be written.
     */
    void beginWrite(double t) override
    {
        if (!multiFile_.is_open()) {
            startMultiFile_(multiFileName_);
        }

        // make sure that all previous output has been written and no other thread
        // accesses the memory used as the target for the extracted quantities
        taskletRunner_.barrier();
        releaseBuffers_();

        curTime_ = t;
        curOutFileName_ = fileName_();

        curWriter_ = std::make_unique<VtkWriter>(gridView_, Dune::VTK::conforming);
        ++curWriterNum_;
    }

    /*!
     * \brief Allocate a managed buffer for a scalar field
     *
     * The buffer will be deleted automatically after the data has
     * been written by to disk.
     */
    ScalarBuffer* allocateManagedScalarBuffer(std::size_t numEntities)
    {
        managedScalarBuffers_.emplace_back(std::make_unique<ScalarBuffer>(numEntities));
        return managedScalarBuffers_.back().get();
    }

    /*!
     * \brief Allocate a managed buffer for a vector field
     *
     * The buffer will be deleted automatically after the data has
     * been written by to disk.
     */
    VectorBuffer* allocateManagedVectorBuffer(std::size_t numOuter, std::size_t numInner)
    {
        managedVectorBuffers_.emplace_back(std::make_unique<VectorBuffer>(numOuter));
        for (auto& buffer : *managedVectorBuffers_.back()) {
            buffer.resize(numInner);
        }

        return managedVectorBuffers_.back().get();
    }

    /*!
     * \brief Add a finished vertex centered vector field to the
     *        output.
     *
     * If the buffer is managed by the VtkMultiWriter, it must have
     * been created using allocateManagedBuffer() and may not be used
     * anywhere after calling this method. After the data is written
     * to disk, it will be deleted automatically.
     *
     * If the buffer is not managed by the MultiWriter, the buffer
     * must exist at least until the call to endWrite()
     * finishes.
     *
     * In both cases, modifying the buffer between the call to this
     * method and endWrite() results in _undefined behavior_.
     */
    void attachScalarVertexData(ScalarBuffer& buf, std::string_view name) override
    {
        sanitizeScalarBuffer_(buf);

        using VtkFn = VtkScalarFunction<GridView, VertexMapper>;
        curWriter_->addVertexData(std::make_shared<VtkFn>(name,
                                                          gridView_,
                                                          vertexMapper_,
                                                          buf,
                                                          /*codim=*/dim));
    }

    /*!
     * \brief Add a element centered quantity to the output.
     *
     * If the buffer is managed by the VtkMultiWriter, it must have
     * been created using createField() and may not be used by
     * anywhere after calling this method. After the data is written
     * to disk, it will be deleted automatically.
     *
     * If the buffer is not managed by the MultiWriter, the buffer
     * must exist at least until the call to endWrite()
     * finishes.
     *
     * In both cases, modifying the buffer between the call to this
     * method and endWrite() results in _undefined behaviour_.
     */
    void attachScalarElementData(ScalarBuffer& buf, std::string_view name) override
    {
        sanitizeScalarBuffer_(buf);

        using VtkFn = VtkScalarFunction<GridView, ElementMapper>;
        curWriter_->addCellData(std::make_shared<VtkFn>(name,
                                                        gridView_,
                                                        elementMapper_,
                                                        buf,
                                                        /*codim=*/0));
    }

    /*!
     * \brief Add a finished vertex centered vector field to the
     *        output.
     *
     * If the buffer is managed by the VtkMultiWriter, it must have
     * been created using allocateManagedBuffer() and may not be used
     * anywhere after calling this method. After the data is written
     * to disk, it will be deleted automatically.
     *
     * If the buffer is not managed by the MultiWriter, the buffer
     * must exist at least until the call to endWrite()
     * finishes.
     *
     * In both cases, modifying the buffer between the call to this
     * method and endWrite() results in _undefined behavior_.
     */
    void attachVectorVertexData(VectorBuffer& buf, std::string_view name) override
    {
        sanitizeVectorBuffer_(buf);

        using VtkFn = VtkVectorFunction<GridView, VertexMapper>;
        curWriter_->addVertexData(std::make_shared<VtkFn>(name,
                                                          gridView_,
                                                          vertexMapper_,
                                                          buf,
                                                          /*codim=*/dim));
    }

    /*!
     * \brief Add a finished vertex-centered tensor field to the output.
     */
    void attachTensorVertexData(TensorBuffer& buf, std::string_view name) override
    {
        using VtkFn = VtkTensorFunction<GridView, VertexMapper>;

        for (unsigned colIdx = 0; colIdx < buf[0].N(); ++colIdx) {
            std::ostringstream oss;
            oss << name <<  "[" << colIdx << "]";
            curWriter_->addVertexData(std::make_shared<VtkFn>(oss.str(),
                                                              gridView_,
                                                              vertexMapper_,
                                                              buf,
                                                              /*codim=*/dim,
                                                              colIdx));
        }
    }

    /*!
     * \brief Add a element centered quantity to the output.
     *
     * If the buffer is managed by the VtkMultiWriter, it must have
     * been created using createField() and may not be used by
     * anywhere after calling this method. After the data is written
     * to disk, it will be deleted automatically.
     *
     * If the buffer is not managed by the MultiWriter, the buffer
     * must exist at least until the call to endWrite()
     * finishes.
     *
     * In both cases, modifying the buffer between the call to this
     * method and endWrite() results in _undefined behaviour_.
     */
    void attachVectorElementData(VectorBuffer& buf, std::string_view name) override
    {
        sanitizeVectorBuffer_(buf);

        using VtkFn = VtkVectorFunction<GridView, ElementMapper>;
        curWriter_->addCellData(std::make_shared<VtkFn>(name,
                                                        gridView_,
                                                        elementMapper_,
                                                        buf,
                                                        /*codim=*/0));
    }

    /*!
     * \brief Add a finished element-centered tensor field to the output.
     */
    void attachTensorElementData(TensorBuffer& buf, std::string_view name) override
    {
        using VtkFn = VtkTensorFunction<GridView, ElementMapper>;

        for (unsigned colIdx = 0; colIdx < buf[0].N(); ++colIdx) {
            std::ostringstream oss;
            oss << name <<  "[" << colIdx << "]";
            curWriter_->addCellData(std::make_shared<VtkFn>(oss.str(),
                                                            gridView_,
                                                            elementMapper_,
                                                            buf,
                                                            /*codim=*/0,
                                                            colIdx));
        }
    }

    /*!
     * \brief Finalizes the current writer.
     *
     * This means that everything will be written to disk, except if
     * the onlyDiscard argument is true. In this case only all managed
     * buffers are deleted, but no output is written.
     */
    void endWrite(bool onlyDiscard) override
    {
        if (!onlyDiscard) {
            auto tasklet = std::make_shared<WriteDataTasklet>(*this);
            taskletRunner_.dispatch(tasklet);
        }
        else {
            --curWriterNum_;
        }

        // temporarily write the closing XML mumbo-jumbo to the mashup
        // file so that the data set can be loaded even if the
        // simulation is aborted (or not yet finished)
        finishMultiFile_();
    }

    /*!
     * \brief Write the multi-writer's state to a restart file.
     */
    template <class Restarter>
    void serialize(Restarter& res)
    {
        res.serializeSectionBegin("VTKMultiWriter");
        res.serializeStream() << curWriterNum_ << "\n";

        if (commRank_ == 0) {
            std::streamsize fileLen = 0;
            std::streamoff filePos = 0;
            if (multiFile_.is_open()) {
                // write the meta file into the restart file
                filePos = multiFile_.tellp();
                multiFile_.seekp(0, std::ios::end);
                fileLen = multiFile_.tellp();
                multiFile_.seekp(filePos);
            }

            res.serializeStream() << fileLen << "  " << filePos << "\n";

            if (fileLen > 0) {
                std::ifstream multiFileIn(multiFileName_.c_str());
                std::vector<char> tmp(fileLen);
                multiFileIn.read(tmp.data(), static_cast<long>(fileLen));
                res.serializeStream().write(tmp.data(), fileLen);
            }
        }

        res.serializeSectionEnd();
    }

    /*!
     * \brief Read the multi-writer's state from a restart file.
     */
    template <class Restarter>
    void deserialize(Restarter& res)
    {
        res.deserializeSectionBegin("VTKMultiWriter");
        res.deserializeStream() >> curWriterNum_;

        if (commRank_ == 0) {
            std::string dummy;
            std::getline(res.deserializeStream(), dummy);

            // recreate the meta file from the restart file
            std::streamoff filePos;
            std::streamsize fileLen;
            res.deserializeStream() >> fileLen >> filePos;
            std::getline(res.deserializeStream(), dummy);
            if (multiFile_.is_open()) {
                multiFile_.close();
            }

            if (fileLen > 0) {
                multiFile_.open(multiFileName_);

                std::vector<char> tmp(fileLen);
                res.deserializeStream().read(tmp.data(), fileLen);
                multiFile_.write(tmp.data(), fileLen);
            }

            multiFile_.seekp(filePos);
        }
        else {
            std::string tmp;
            std::getline(res.deserializeStream(), tmp);
        }
        res.deserializeSectionEnd();
    }

private:
    std::string fileName_() const
    {
        // use a new file name for each time step
        std::ostringstream oss;
        oss << simName_ << "-"
            << std::setw(5) << std::setfill('0') << curWriterNum_;
        return oss.str();
    }

    static std::string fileSuffix_()
    { return (GridView::dimension == 1) ? "vtp" : "vtu"; }

    void startMultiFile_(const std::string& multiFileName)
    {
        // only the first process writes to the multi-file
        if (commRank_ == 0) {
            // generate one meta vtk-file holding the individual time steps
            multiFile_.open(multiFileName);
            multiFile_ << "<?xml version=\"1.0\"?>\n"
                          "<VTKFile type=\"Collection\"\n"
                          "         version=\"0.1\"\n"
                          "         byte_order=\"LittleEndian\"\n"
                          "         compressor=\"vtkZLibDataCompressor\">\n"
                          " <Collection>\n";
        }
    }

    void finishMultiFile_()
    {
        // only the first process writes to the multi-file
        if (commRank_ == 0) {
            // make sure that we always have a working meta file
            std::ofstream::pos_type pos = multiFile_.tellp();
            multiFile_ << " </Collection>\n"
                          "</VTKFile>\n";
            multiFile_.seekp(pos);
            multiFile_.flush();
        }
    }

    // make sure the field is well defined if running under valgrind
    // and make sure that all values can be displayed by paraview
    void sanitizeScalarBuffer_(ScalarBuffer&)
    {
        // nothing to do: this is done by VtkScalarFunction
    }

    void sanitizeVectorBuffer_(VectorBuffer&)
    {
        // nothing to do: this is done by VtkVectorFunction
    }

    // release the memory occupied by all buffer objects managed by the multi-writer
    void releaseBuffers_()
    {
        // discard managed objects and the current VTK writer
        curWriter_.reset();
        managedScalarBuffers_.clear();
        managedVectorBuffers_.clear();
    }

    const GridView gridView_;
    ElementMapper elementMapper_;
    VertexMapper vertexMapper_;

    const std::string outputDir_;
    const std::string simName_;
    std::ofstream multiFile_;
    const std::string multiFileName_;

    const int commSize_; // number of processes in the communicator
    const int commRank_; // rank of the current process in the communicator

    std::unique_ptr<VtkWriter> curWriter_;
    double curTime_{};
    std::string curOutFileName_;
    int curWriterNum_;

    std::vector<std::unique_ptr<ScalarBuffer>> managedScalarBuffers_;
    std::vector<std::unique_ptr<VectorBuffer>> managedVectorBuffers_;

    TaskletRunner taskletRunner_;
};

} // namespace Opm

#endif
