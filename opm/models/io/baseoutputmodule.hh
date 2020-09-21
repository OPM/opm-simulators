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
 * \copydoc Opm::BaseOutputModule
 */
#ifndef EWOMS_BASE_OUTPUT_MODULE_HH
#define EWOMS_BASE_OUTPUT_MODULE_HH

#include "baseoutputwriter.hh"

#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/basicproperties.hh>
#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>

#include <vector>
#include <sstream>
#include <string>
#include <array>

#include <cstdio>

namespace Opm::Properties {

template <class TypeTag, class MyTypeTag>
struct FluidSystem;

} // namespace Opm::Properties

namespace Opm {

#if __GNUC__ || __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
#endif

/*!
 * \brief The base class for writer modules.
 *
 * This class also provides some convenience methods for buffer
 * management and is the base class for all other output writer
 * modules.
 */
template<class TypeTag>
class BaseOutputModule
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using DiscBaseOutputModule = GetPropType<TypeTag, Properties::DiscBaseOutputModule>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    using Tensor = BaseOutputWriter::Tensor;

public:
    using ScalarBuffer = BaseOutputWriter::ScalarBuffer;
    using VectorBuffer = BaseOutputWriter::VectorBuffer;
    using TensorBuffer = BaseOutputWriter::TensorBuffer;

    using EqBuffer = std::array<ScalarBuffer, numEq>;
    using PhaseBuffer = std::array<ScalarBuffer, numPhases>;
    using ComponentBuffer = std::array<ScalarBuffer, numComponents>;
    using PhaseComponentBuffer = std::array<std::array<ScalarBuffer, numComponents>, numPhases>;

    using PhaseVectorBuffer = std::array<VectorBuffer, numPhases>;

    BaseOutputModule(const Simulator& simulator)
        : simulator_(simulator)
    {}

    virtual ~BaseOutputModule()
    {}

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to disk.
     *
     * The module can dynamically cast the writer to the desired
     * concrete class. If the writer is incompatible with the module,
     * this method should become a no-op.
     */
    virtual void allocBuffers() = 0;

    /*!
     * \brief Modify the internal buffers according to the intensive quanties relevant
     *        for an element
     *
     * The module can dynamically cast the writer to the desired
     * concrete class. If the writer is incompatible with the module,
     * this method should become a no-op.
     */
    virtual void processElement(const ElementContext& elemCtx) = 0;

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    virtual void commitBuffers(BaseOutputWriter& writer) = 0;

    /*!
     * \brief Returns true iff the module needs to access the extensive quantities of a
     * context to do its job.
     *
     * For example, this happens if velocities or gradients should be written.
     *
     * Always returning true here does not do any harm from the correctness perspective,
     * but it slows down writing the output fields. Since most output modules only write
     * intensive quantities, this method returns 'false' by default.
     */
    virtual bool needExtensiveQuantities() const
    { return false; }

protected:
    enum BufferType {
        //! Buffer contains data associated with the degrees of freedom
        DofBuffer,

        //! Buffer contains data associated with the grid's vertices
        VertexBuffer,

        //! Buffer contains data associated with the grid's elements
        ElementBuffer
    };

    /*!
     * \brief Allocate the space for a buffer storing a scalar quantity
     */
    void resizeScalarBuffer_(ScalarBuffer& buffer,
                             BufferType bufferType = DofBuffer)
    {
        size_t n;
        if (bufferType == VertexBuffer)
            n = static_cast<size_t>(simulator_.gridView().size(dim));
        else if (bufferType == ElementBuffer)
            n = static_cast<size_t>(simulator_.gridView().size(0));
        else if (bufferType == DofBuffer)
            n = simulator_.model().numGridDof();
        else
            throw std::logic_error("bufferType must be one of Dof, Vertex or Element");

        buffer.resize(n);
        std::fill(buffer.begin(), buffer.end(), 0.0);
    }

    /*!
     * \brief Allocate the space for a buffer storing a tensorial quantity
     */
    void resizeTensorBuffer_(TensorBuffer& buffer,
                             BufferType bufferType = DofBuffer)
    {
        size_t n;
        if (bufferType == VertexBuffer)
            n = static_cast<size_t>(simulator_.gridView().size(dim));
        else if (bufferType == ElementBuffer)
            n = static_cast<size_t>(simulator_.gridView().size(0));
        else if (bufferType == DofBuffer)
            n = simulator_.model().numGridDof();
        else
            throw std::logic_error("bufferType must be one of Dof, Vertex or Element");

        buffer.resize(n);
        Tensor nullMatrix(dimWorld, dimWorld, 0.0);
        std::fill(buffer.begin(), buffer.end(), nullMatrix);
    }

    /*!
     * \brief Allocate the space for a buffer storing a equation specific
     *        quantity
     */
    void resizeEqBuffer_(EqBuffer& buffer,
                         BufferType bufferType = DofBuffer)
    {
        size_t n;
        if (bufferType == VertexBuffer)
            n = static_cast<size_t>(simulator_.gridView().size(dim));
        else if (bufferType == ElementBuffer)
            n = static_cast<size_t>(simulator_.gridView().size(0));
        else if (bufferType == DofBuffer)
            n = simulator_.model().numGridDof();
        else
            throw std::logic_error("bufferType must be one of Dof, Vertex or Element");

        for (unsigned i = 0; i < numEq; ++i) {
            buffer[i].resize(n);
            std::fill(buffer[i].begin(), buffer[i].end(), 0.0);
        }
    }

    /*!
     * \brief Allocate the space for a buffer storing a phase-specific
     *        quantity
     */
    void resizePhaseBuffer_(PhaseBuffer& buffer,
                            BufferType bufferType = DofBuffer)
    {
        size_t n;
        if (bufferType == VertexBuffer)
            n = static_cast<size_t>(simulator_.gridView().size(dim));
        else if (bufferType == ElementBuffer)
            n = static_cast<size_t>(simulator_.gridView().size(0));
        else if (bufferType == DofBuffer)
            n = simulator_.model().numGridDof();
        else
            throw std::logic_error("bufferType must be one of Dof, Vertex or Element");

        for (unsigned i = 0; i < numPhases; ++i) {
            buffer[i].resize(n);
            std::fill(buffer[i].begin(), buffer[i].end(), 0.0);
        }
    }

    /*!
     * \brief Allocate the space for a buffer storing a component
     *        specific quantity
     */
    void resizeComponentBuffer_(ComponentBuffer& buffer,
                                BufferType bufferType = DofBuffer)
    {
        size_t n;
        if (bufferType == VertexBuffer)
            n = static_cast<size_t>(simulator_.gridView().size(dim));
        else if (bufferType == ElementBuffer)
            n = static_cast<size_t>(simulator_.gridView().size(0));
        else if (bufferType == DofBuffer)
            n = simulator_.model().numGridDof();
        else
            throw std::logic_error("bufferType must be one of Dof, Vertex or Element");

        for (unsigned i = 0; i < numComponents; ++i) {
            buffer[i].resize(n);
            std::fill(buffer[i].begin(), buffer[i].end(), 0.0);
        }
    }

    /*!
     * \brief Allocate the space for a buffer storing a phase and
     *        component specific buffer
     */
    void resizePhaseComponentBuffer_(PhaseComponentBuffer& buffer,
                                     BufferType bufferType = DofBuffer)
    {
        size_t n;
        if (bufferType == VertexBuffer)
            n = static_cast<size_t>(simulator_.gridView().size(dim));
        else if (bufferType == ElementBuffer)
            n = static_cast<size_t>(simulator_.gridView().size(0));
        else if (bufferType == DofBuffer)
            n = simulator_.model().numGridDof();
        else
            throw std::logic_error("bufferType must be one of Dof, Vertex or Element");

        for (unsigned i = 0; i < numPhases; ++i) {
            for (unsigned j = 0; j < numComponents; ++j) {
                buffer[i][j].resize(n);
                std::fill(buffer[i][j].begin(), buffer[i][j].end(), 0.0);
            }
        }
    }

    /*!
     * \brief Add a buffer containing scalar quantities to the result file.
     */
    void commitScalarBuffer_(BaseOutputWriter& baseWriter,
                             const char *name,
                             ScalarBuffer& buffer,
                             BufferType bufferType = DofBuffer)
    {
        if (bufferType == DofBuffer)
            DiscBaseOutputModule::attachScalarDofData_(baseWriter, buffer, name);
        else if (bufferType == VertexBuffer)
            attachScalarVertexData_(baseWriter, buffer, name);
        else if (bufferType == ElementBuffer)
            attachScalarElementData_(baseWriter, buffer, name);
        else
            throw std::logic_error("bufferType must be one of Dof, Vertex or Element");
    }

    /*!
     * \brief Add a buffer containing vectorial quantities to the result file.
     */
    void commitVectorBuffer_(BaseOutputWriter& baseWriter,
                             const char *name,
                             VectorBuffer& buffer,
                             BufferType bufferType = DofBuffer)
    {
        if (bufferType == DofBuffer)
            DiscBaseOutputModule::attachVectorDofData_(baseWriter, buffer, name);
        else if (bufferType == VertexBuffer)
            attachVectorVertexData_(baseWriter, buffer, name);
        else if (bufferType == ElementBuffer)
            attachVectorElementData_(baseWriter, buffer, name);
        else
            throw std::logic_error("bufferType must be one of Dof, Vertex or Element");
    }

    /*!
     * \brief Add a buffer containing tensorial quantities to the result file.
     */
    void commitTensorBuffer_(BaseOutputWriter& baseWriter,
                             const char *name,
                             TensorBuffer& buffer,
                             BufferType bufferType = DofBuffer)
    {
        if (bufferType == DofBuffer)
            DiscBaseOutputModule::attachTensorDofData_(baseWriter, buffer, name);
        else if (bufferType == VertexBuffer)
            attachTensorVertexData_(baseWriter, buffer, name);
        else if (bufferType == ElementBuffer)
            attachTensorElementData_(baseWriter, buffer, name);
        else
            throw std::logic_error("bufferType must be one of Dof, Vertex or Element");
    }

    /*!
     * \brief Add a buffer with as many variables as PDEs to the result file.
     */
    void commitPriVarsBuffer_(BaseOutputWriter& baseWriter,
                              const char *pattern,
                              EqBuffer& buffer,
                              BufferType bufferType = DofBuffer)
    {
        char name[512];
        for (unsigned i = 0; i < numEq; ++i) {
            std::string eqName = simulator_.model().primaryVarName(i);
            snprintf(name, 512, pattern, eqName.c_str());

            if (bufferType == DofBuffer)
                DiscBaseOutputModule::attachScalarDofData_(baseWriter, buffer[i], name);
            else if (bufferType == VertexBuffer)
                attachScalarVertexData_(baseWriter, buffer[i], name);
            else if (bufferType == ElementBuffer)
                attachScalarElementData_(baseWriter, buffer[i], name);
            else
                throw std::logic_error("bufferType must be one of Dof, Vertex or Element");
        }
    }

    /*!
     * \brief Add a buffer with as many variables as PDEs to the result file.
     */
    void commitEqBuffer_(BaseOutputWriter& baseWriter,
                         const char *pattern,
                         EqBuffer& buffer,
                         BufferType bufferType = DofBuffer)
    {
        char name[512];
        for (unsigned i = 0; i < numEq; ++i) {
            std::ostringstream oss;
            oss << i;
            snprintf(name, 512, pattern, oss.str().c_str());

            if (bufferType == DofBuffer)
                DiscBaseOutputModule::attachScalarDofData_(baseWriter, buffer[i], name);
            else if (bufferType == VertexBuffer)
                attachScalarVertexData_(baseWriter, buffer[i], name);
            else if (bufferType == ElementBuffer)
                attachScalarElementData_(baseWriter, buffer[i], name);
            else
                throw std::logic_error("bufferType must be one of Dof, Vertex or Element");
        }
    }

    /*!
     * \brief Add a phase-specific buffer to the result file.
     */
    void commitPhaseBuffer_(BaseOutputWriter& baseWriter,
                            const char *pattern,
                            PhaseBuffer& buffer,
                            BufferType bufferType = DofBuffer)
    {
        char name[512];
        for (unsigned i = 0; i < numPhases; ++i) {
            snprintf(name, 512, pattern, FluidSystem::phaseName(i));

            if (bufferType == DofBuffer)
                DiscBaseOutputModule::attachScalarDofData_(baseWriter, buffer[i], name);
            else if (bufferType == VertexBuffer)
                attachScalarVertexData_(baseWriter, buffer[i], name);
            else if (bufferType == ElementBuffer)
                attachScalarElementData_(baseWriter, buffer[i], name);
            else
                throw std::logic_error("bufferType must be one of Dof, Vertex or Element");
        }
    }

    /*!
     * \brief Add a component-specific buffer to the result file.
     */
    void commitComponentBuffer_(BaseOutputWriter& baseWriter,
                                const char *pattern,
                                ComponentBuffer& buffer,
                                BufferType bufferType = DofBuffer)
    {
        char name[512];
        for (unsigned i = 0; i < numComponents; ++i) {
            snprintf(name, 512, pattern, FluidSystem::componentName(i));

            if (bufferType == DofBuffer)
                DiscBaseOutputModule::attachScalarDofData_(baseWriter, buffer[i], name);
            else if (bufferType == VertexBuffer)
                attachScalarVertexData_(baseWriter, buffer[i], name);
            else if (bufferType == ElementBuffer)
                attachScalarElementData_(baseWriter, buffer[i], name);
            else
                throw std::logic_error("bufferType must be one of Dof, Vertex or Element");
        }
    }

    /*!
     * \brief Add a phase and component specific quantities to the output.
     */
    void commitPhaseComponentBuffer_(BaseOutputWriter& baseWriter,
                                     const char *pattern,
                                     PhaseComponentBuffer& buffer,
                                     BufferType bufferType = DofBuffer)
    {
        char name[512];
        for (unsigned i= 0; i < numPhases; ++i) {
            for (unsigned j = 0; j < numComponents; ++j) {
                snprintf(name, 512, pattern,
                         FluidSystem::phaseName(i),
                         FluidSystem::componentName(j));

                if (bufferType == DofBuffer)
                    DiscBaseOutputModule::attachScalarDofData_(baseWriter, buffer[i][j], name);
                else if (bufferType == VertexBuffer)
                    attachScalarVertexData_(baseWriter, buffer[i][j], name);
                else if (bufferType == ElementBuffer)
                    attachScalarElementData_(baseWriter, buffer[i][j], name);
                else
                    throw std::logic_error("bufferType must be one of Dof, Vertex or Element");
            }
        }
    }

    void attachScalarElementData_(BaseOutputWriter& baseWriter,
                                  ScalarBuffer& buffer,
                                  const char *name)
    { baseWriter.attachScalarElementData(buffer, name); }

    void attachScalarVertexData_(BaseOutputWriter& baseWriter,
                                 ScalarBuffer& buffer,
                                 const char *name)
    { baseWriter.attachScalarVertexData(buffer, name); }

    void attachVectorElementData_(BaseOutputWriter& baseWriter,
                                  VectorBuffer& buffer,
                                  const char *name)
    { baseWriter.attachVectorElementData(buffer, name); }

    void attachVectorVertexData_(BaseOutputWriter& baseWriter,
                                 VectorBuffer& buffer,
                                 const char *name)
    { baseWriter.attachVectorVertexData(buffer, name); }

    void attachTensorElementData_(BaseOutputWriter& baseWriter,
                                  TensorBuffer& buffer,
                                  const char *name)
    { baseWriter.attachTensorElementData(buffer, name); }

    void attachTensorVertexData_(BaseOutputWriter& baseWriter,
                                 TensorBuffer& buffer,
                                 const char *name)
    { baseWriter.attachTensorVertexData(buffer, name); }

    const Simulator& simulator_;
};

#if __GNUC__ || __clang__
#pragma GCC diagnostic pop
#endif

} // namespace Opm

#endif
