# Porting a Composite Class to the GPU

This document describes how to port an existing composite C++ class
(i.e. a class that *owns* one or more containers such as `std::vector`)
to the GPU using the **`GpuBuffer`** / **`GpuView`** paradigm that is
already established in `opm-simulators/opm/simulators/linalg/gpuistl/`.

> **Golden rule: edit, do not duplicate.**
> The goal is to add GPU support *to the existing class* by templating
> it on its storage container. **Do not create a separate GPU class.**
> A duplicate class doubles the maintenance burden, risks divergence,
> and is rejected in code review. Creating a brand-new GPU-only class
> is an exceptional last resort reserved for situations where the
> original class truly cannot be templated (e.g. it has a heavily
> virtual interface or is generated code). In all normal cases the
> single header file is edited in-place.

The pattern is split across three repositories:

| Repository      | Role                                                                                  |
| --------------- | ------------------------------------------------------------------------------------- |
| `opm-common`    | Header-only physics / material classes that need to run on the GPU.                   |
| `opm-grid`      | Header-only grid / utility classes that need to run on the GPU (e.g. `SparseTable`).  |
| `opm-simulators`| Hosts `GpuBuffer`, `GpuView`, the CUDA/HIP build system, the kernels and **all** tests.|

> **Important:** Only `opm-simulators` actually compiles GPU code.
> `opm-common` and `opm-grid` may only contain *header-only* GPU code,
> and any GPU-specific helper (such as `copy_to_gpu` / `make_view`)
> **must** be guarded by `#if HAVE_CUDA` (see `SparseTable.hpp` for an
> example). The class itself stays GPU-agnostic; the GPU instantiation
> happens through a template parameter for the storage container.

---

## 1. The two pillars: `GpuBuffer` and `GpuView`

`opm-simulators/opm/simulators/linalg/gpuistl/`

* **`GpuBuffer<T>`** — *owning* container. It allocates GPU memory in
  its constructor and frees it in its destructor (RAII). It is the
  GPU equivalent of `std::vector<T>` and lives on the **host** (the
  object itself is a host object, but the bytes it points to are on
  the device).
* **`GpuView<T>`** — *non-owning* view over a contiguous block of GPU
  memory. It is essentially `{T* ptr; size_t size;}` and is **trivially
  copyable**, which is what makes it safe to pass *by value* into a
  `__global__` kernel.

The relationship is intentionally analogous to
`std::vector<T>` ↔ `std::span<T>`.

### Ownership rules (read this before writing any code)

1. **`copy_to_gpu` returns an owning object built from `GpuBuffer`s.**
   The caller (CPU code) **must keep this object alive** for as long as
   the GPU is allowed to dereference the data. Destroying it frees the
   GPU memory.
2. **`make_view` returns a non-owning object built from `GpuView`s.**
   The view object is the thing you pass into kernels. It does **not**
   extend the lifetime of the underlying buffers, so the corresponding
   buffer-holding object from step 1 must outlive every kernel that
   uses the view.

A typical CPU-side flow looks like this:

```cpp
// 1. CPU object, fully populated and finalized.
auto cpuParams = makeMyParams();

// 2. Owning GPU copy — return type is always MyParams<..., GpuBuffer>, no template arg needed.
//    Hold on to this for the entire lifetime of all kernel launches that need the data.
auto gpuBufferParams = Opm::gpuistl::copy_to_gpu(cpuParams);

// 3. Non-owning view — return type is always MyParams<..., GpuView>, no template arg needed.
//    Cheap to copy, safe to pass to kernels.
auto gpuViewParams = Opm::gpuistl::make_view(gpuBufferParams);

myKernel<<<grid, block>>>(gpuViewParams, /* … */);   // pass by value
OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());

// 4. gpuBufferParams goes out of scope here → device memory is freed.
```

---

## 2. Decorating the class: `gpuDecorators.hpp`

`opm-common/opm/common/utility/gpuDecorators.hpp`

> When creating a GPU-compatible class, **use `OPM_HOST_DEVICE` instead
> of `__host__ __device__`**. Include the header:
>
> ```cpp
> #include <opm/common/utility/gpuDecorators.hpp>
> ```

The header is safe to include unconditionally — when the GPU compiler
is not active, the macros expand to nothing. The full list of macros it
provides is:

| Macro                                | Expands to (GPU build)        | Expands to (CPU build) | Purpose                                                           |
| ------------------------------------ | ----------------------------- | ---------------------- | ----------------------------------------------------------------- |
| `OPM_HOST_DEVICE`                    | `__device__ __host__`         | *(empty)*              | Decorate functions callable from both CPU and a CUDA/HIP kernel.  |
| `OPM_DEVICE`                         | `__device__`                  | *(empty)*              | Decorate functions only callable from the device.                 |
| `OPM_HOST`                           | `__host__`                    | *(empty)*              | Included for completeness; **rarely needed** in practice — most host-only logic is compiled by the host compiler without any decorator. |
| `OPM_IS_USING_GPU`                   | `1`                           | `0`                    | Compile-time constant: GPU support enabled (`HAVE_CUDA`) or not.  |
| `OPM_IS_COMPILING_WITH_GPU_COMPILER` | `1` under `nvcc` / `hipcc`    | `0`                    | True when the current TU is being processed by `nvcc` or `hipcc`. |
| `OPM_IS_INSIDE_DEVICE_FUNCTION`      | `1` inside a device function  | `0`                    | Useful for `if constexpr`-style specialization within a function. |
| `OPM_IS_INSIDE_HOST_FUNCTION`        | `0` inside a device function  | `1`                    | The complement of the above.                                      |

When `HAVE_CUDA` is set the header also pulls in either
`<cuda_runtime.h>` or `<hip/hip_runtime.h>` (selected via `USE_HIP`).

### Exceptions: use `OPM_THROW`

Always use **`OPM_THROW`** (from `<opm/common/ErrorMacros.hpp>`) for
error reporting; do not use raw `throw`. `OPM_THROW` is GPU-aware and
behaves correctly in code shared between host and device:

```cpp
#include <opm/common/ErrorMacros.hpp>

OPM_HOST_DEVICE double safeDiv(double a, double b)
{
    if (b == 0.0) {
        OPM_THROW(std::runtime_error, "Division by zero");
    }
    return a / b;
}
```

There is a dedicated test (`tests/gpuistl/test_throw_macros_on_gpu.cu`)
that validates this behaviour.

---

## 3. Templating the class on its storage container

> **Edit the existing file — do not copy it.**
> Open the class's existing `.hpp` file and add the `Storage` template
> parameter directly. The CPU instantiation (`Storage =
> VectorWithDefaultAllocator`) must continue to behave exactly as
> before — callers that do not supply the second template argument are
> completely unaffected.

The trick that makes the same class instantiable on CPU *and* GPU is
to template it on the container type using a **template template
parameter**, defaulting to `VectorWithDefaultAllocator`.
The same class is then instantiated three times:

| Instantiation                        | Storage type                  | Lives on |
| ------------------------------------ | ----------------------------- | -------- |
| `MyClass<Traits>`                    | `VectorWithDefaultAllocator<Scalar>` | CPU |
| `MyClass<Traits, GpuBuffer>`         | `GpuBuffer<Scalar>` (owning)  | GPU mem, CPU object |
| `MyClass<Traits, GpuView>`           | `GpuView<Scalar>`             | GPU (passed to kernel) |

**Limit the public interface for GPU instantiations.**
When `Storage` is `GpuBuffer` or `GpuView`, many methods that are
perfectly valid on the CPU (e.g. those that resize, append, serialize,
or return `std::string`) are meaningless or impossible on the device.
Do not try to make every method work for every storage type. Instead:

* Guard host-only methods with `if constexpr` or a `static_assert` so
  that they produce a clear compile-time error when called on a
  `GpuBuffer` / `GpuView` instantiation.
* Decorate only the methods that *need* to run in a kernel with
  `OPM_HOST_DEVICE`; leave the rest undecorated.
* It is acceptable — and encouraged — for `MyClass<Traits, GpuView>`
  to expose fewer methods than `MyClass<Traits>`. The GPU instantiation
  is a read-only, compute-only slice of the full class.

The canonical reference implementation is
`opm-common/opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterialParams.hpp`.
Other examples worth studying:

* `opm-common/opm/material/common/UniformTabulated2DFunction.hpp`
* `opm-common/opm/material/components/CO2Tables.hpp`
* `opm-common/opm/material/fluidsystems/blackoilpvt/BrineCo2Pvt.hpp`
* `opm-grid/opm/grid/utility/SparseTable.hpp`

### Skeleton of a GPU-portable composite class

```cpp
// my_composite.hpp  (lives in opm-common or opm-grid; HEADER ONLY)
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/common/utility/VectorWithDefaultAllocator.hpp>

// Use template<class> class Storage (template template parameter) so that
// GpuBuffer, GpuView, and VectorWithDefaultAllocator can all be passed
// without spelling out the scalar type at the call site.
template <class Traits, template<class> class Storage = VectorWithDefaultAllocator>
class MyComposite;

// Forward-declare the GPU helpers in their final namespace so that we
// can friend them inside the class.
// Both are guarded by HAVE_CUDA; they always return GpuBuffer / GpuView
// instantiations — the caller never needs to spell out the storage type.
namespace Opm::gpuistl {
    template <class Traits>
    MyComposite<Traits, GpuBuffer>
    copy_to_gpu(const MyComposite<Traits>& cpu);

    template <class Traits, template<class> class ContainerT>
    MyComposite<Traits, GpuView>
    make_view(MyComposite<Traits, ContainerT>& gpuBuffers);
}

namespace Opm {

template <class Traits, template<class> class Storage>
class MyComposite
{
    using Scalar  = typename Traits::Scalar;
    using VectorT = Storage<Scalar>;   // concrete type for this instantiation
public:
    MyComposite() = default;
    MyComposite(VectorT xs, VectorT ys)
        : xs_(std::move(xs)), ys_(std::move(ys))
    {
        if (xs_.size() != ys_.size()) {
            OPM_THROW(std::invalid_argument,
                      "xs and ys must have equal size");
        }
    }

    OPM_HOST_DEVICE const VectorT& xs() const { return xs_; }
    OPM_HOST_DEVICE const VectorT& ys() const { return ys_; }

    OPM_HOST_DEVICE Scalar evaluate(Scalar x) const
    {
        // ... uses xs_/ys_ via operator[] which is also OPM_HOST_DEVICE
    }

private:
    // Friends must match the forward declarations exactly.
    template <class T>
    friend MyComposite<T, GpuBuffer>
    Opm::gpuistl::copy_to_gpu(const MyComposite<T>& cpu);

    template <class T, template<class> class ContainerT>
    friend MyComposite<T, GpuView>
    Opm::gpuistl::make_view(MyComposite<T, ContainerT>& gpuBuffers);

    VectorT xs_;
    VectorT ys_;
};

} // namespace Opm
```

---

## 4. Writing `copy_to_gpu`

`copy_to_gpu` returns an **owning** GPU object. It is the moment when
host memory is transferred to the device.

```cpp
// Still in the same header, but inside namespace Opm::gpuistl and
// guarded by HAVE_CUDA in repos that do not compile GPU code:
#if HAVE_CUDA
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>

namespace Opm::gpuistl {

template <class Traits>
MyComposite<Traits, GpuBuffer>
copy_to_gpu(const MyComposite<Traits>& cpu)
{
    using Scalar = typename Traits::Scalar;

    // Defensive checks on the host object first.
    if (cpu.xs().size() != cpu.ys().size()) {
        OPM_THROW(std::logic_error,
                  "Inconsistent CPU object passed to copy_to_gpu");
    }

    // Return type is always MyComposite<Traits, GpuBuffer>.
    // Each constructor performs a synchronous H→D copy.
    auto gpuXs = GpuBuffer<Scalar>(cpu.xs());
    auto gpuYs = GpuBuffer<Scalar>(cpu.ys());

    return MyComposite<Traits, GpuBuffer>(gpuXs, gpuYs);
}

} // namespace Opm::gpuistl
#endif // HAVE_CUDA
```

Key points:

* **Return type:** `copy_to_gpu` always returns
  `MyClass<OtherTemplates, GpuBuffer>` — the storage type is fixed, not
  a template argument. **Use `auto` to hold the result** and call with
  no explicit template argument:
  ```cpp
  // No template argument needed — return type is always MyClass<Traits, GpuBuffer>
  auto owned = Opm::gpuistl::copy_to_gpu(cpuObj);
  ```
* The returned object **owns** every byte of GPU memory it points to,
  via the `GpuBuffer` members. Keep it alive on the CPU for as long as
  any kernel may dereference it.
* Use `OPM_THROW` for any precondition violation. Heavyweight checks
  (e.g. "is the object finalized?") belong here, not in the kernel.

---

## 5. Writing `make_view`

`make_view` returns a **non-owning** object. It exists solely to be
passed *by value* into a kernel.

```cpp
#if HAVE_CUDA
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>

namespace Opm::gpuistl {

template <class Traits, template<class> class ContainerT>
MyComposite<Traits, GpuView>
make_view(MyComposite<Traits, ContainerT>& gpuBuffers)
{
    using Scalar = typename Traits::Scalar;

    // Return type is always MyComposite<Traits, GpuView>.
    // make_view(GpuBuffer<T>&) is provided by GpuBuffer.hpp.
    GpuView<Scalar> xsView = make_view<Scalar>(gpuBuffers.xs_);
    GpuView<Scalar> ysView = make_view<Scalar>(gpuBuffers.ys_);

    return MyComposite<Traits, GpuView>(xsView, ysView);
}

} // namespace Opm::gpuistl
#endif // HAVE_CUDA
```

Key points:

* **Return type:** `make_view` always returns
  `MyClass<OtherTemplates, GpuView>` — the storage type is fixed, not
  a template argument. This type is trivially copyable and safe to pass
  *by value* to a `__global__` kernel. **Use `auto` to hold the result**
  and call with no explicit template argument:
  ```cpp
  // No template argument needed — return type is always MyClass<Traits, GpuView>
  auto view = Opm::gpuistl::make_view(owned);
  ```
* `ViewT` is typically `GpuView`. The `const`-ness of the element type
  is baked into `GpuView` itself for read-only kernels, which is the
  strongest signal of that intent.
* Because it is non-owning, the original `gpuBuffers` object **must
  outlive every kernel that uses the view**.
* `make_view` accesses private members directly (hence the `friend`
  declarations in section 3). For primitive scalars there is no work
  to do other than wrapping pointers; no H↔D transfer happens here.

---

## 5b. GPU smart pointers (`gpu_smart_pointer.hpp`)

`opm-simulators/opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp`

When a class member is a *scalar or small struct* (rather than a
vector), you may want to store it in GPU memory via a smart pointer
rather than wrapping it in `GpuBuffer`. The utilities in
`gpu_smart_pointer.hpp` handle this.

### Allocation and copying

| Function | Description |
| -------- | ----------- |
| `make_gpu_shared_ptr<T>()` | Allocates `sizeof(T)` bytes on the GPU; returns `std::shared_ptr<T>` with a `cudaFree` deleter. |
| `make_gpu_shared_ptr<T>(value)` | Same, but immediately copies `value` host → device. |
| `make_gpu_unique_ptr<T>()` | Same idea, returns `std::unique_ptr<T, ...>`. |
| `make_gpu_unique_ptr<T>(value)` | Same, with initial host → device copy. |
| `copyFromGPU(ptr)` | Synchronous device → host copy; accepts a raw pointer, `shared_ptr`, or `unique_ptr`. |
| `copyToGPU(value, ptr)` | Synchronous host → device copy; accepts the same pointer types. |

```cpp
// Allocate and initialise a single Scalar on the GPU.
auto gpuVal = Opm::gpuistl::make_gpu_unique_ptr<Scalar>(hostValue);

// Later, read it back (involves a synchronization point).
Scalar result = Opm::gpuistl::copyFromGPU(gpuVal);
```

### Passing smart-pointer data to kernels: `PointerView`

`std::shared_ptr` and `std::unique_ptr` are **not** trivially copyable,
so they cannot be passed by value to a `__global__` kernel. Use
`PointerView<T>` instead — it is a lightweight, trivially copyable
wrapper around a raw `T*`:

```cpp
// make_view(smart_ptr) returns PointerView<T> — non-owning and
// trivially copyable, safe to pass by value to a kernel.
auto view = Opm::gpuistl::make_view(gpuVal);   // PointerView<Scalar>

myKernel<<<grid, block>>>(view, /* … */);       // pass by value
```

Inside the kernel `PointerView<T>` behaves exactly like a pointer:

```cpp
__global__
void myKernel(Opm::gpuistl::PointerView<Scalar> val, Scalar* out)
{
    *out = *val + 1.0f;        // dereference
    // or: *out = val.get()[0];
}
```

The lifetime rule mirrors `GpuView`: the underlying smart pointer
**must outlive** every kernel that uses the view.

### Wrapping an inline value: `ValueAsPointer`

`ValueAsPointer<T>` wraps a plain `T` so that it can be used wherever
a pointer-like interface is expected inside a kernel. This is useful
for generic template code that must work with both `PointerView<T>`
(pointing at GPU memory) and an inline value (no allocation):

```cpp
Opm::gpuistl::ValueAsPointer<float> vp(3.14f);

__global__
void myKernel(Opm::gpuistl::ValueAsPointer<float> v, float* out)
{
    *out = *v;    // dereferences like a pointer; no GPU allocation needed
}
```

---

## 5c. Convenience include: `gpuistl_if_available.hpp`

`opm-common/opm/common/utility/gpuistl_if_available.hpp`

Sections 4 and 5 showed `copy_to_gpu` / `make_view` bodies that manually
include individual GPU headers inside `#if HAVE_CUDA` / `#if USE_HIP`
blocks. In practice, replace all of that boilerplate with a single include
of **`gpuistl_if_available.hpp`**. The header does everything for you:

* It always includes `gpuDecorators.hpp` (so `OPM_HOST_DEVICE` etc. are
  available unconditionally).
* When `HAVE_CUDA` is defined it selects the correct backend:
  * `USE_HIP` → `gpuistl_hip/GpuBuffer.hpp`, `GpuView.hpp`,
    `GpuVector.hpp`, `gpu_smart_pointer.hpp`
  * otherwise → the corresponding `gpuistl/` headers
* When `HAVE_CUDA` is **not** defined the header is empty (beyond the
  `gpuDecorators.hpp` pull-in).

### Usage in `opm-common`

Include unconditionally at the top of the header — no `#if HAVE_CUDA`
guard is required:

```cpp
// In an opm-common header:
#include <opm/common/utility/gpuistl_if_available.hpp>

// copy_to_gpu / make_view still need their own #if HAVE_CUDA guard
// because the function *bodies* reference GpuBuffer/GpuView types,
// but the include itself does not:
#if HAVE_CUDA
namespace Opm::gpuistl {

template <class Traits>
MyComposite<Traits, GpuBuffer>
copy_to_gpu(const MyComposite<Traits>& cpu)
{
    // GpuBuffer is already available — no additional #include needed.
    return MyComposite<Traits, GpuBuffer>(
        GpuBuffer<typename Traits::Scalar>(cpu.xs()),
        GpuBuffer<typename Traits::Scalar>(cpu.ys()));
}

} // namespace Opm::gpuistl
#endif // HAVE_CUDA
```

### Usage in `opm-grid`

`opm-grid` may be built without `opm-common` as a dependency. Guard the
include with `HAVE_OPM_COMMON`:

```cpp
// In an opm-grid header:
#if HAVE_OPM_COMMON
#include <opm/common/utility/gpuistl_if_available.hpp>
#endif
```

The `copy_to_gpu` / `make_view` bodies are then double-guarded:

```cpp
#if HAVE_OPM_COMMON && HAVE_CUDA
namespace Opm::gpuistl {
// ... implementations using GpuBuffer / GpuView ...
} // namespace Opm::gpuistl
#endif // HAVE_OPM_COMMON && HAVE_CUDA
```

See `opm-grid/opm/grid/utility/SparseTable.hpp` for a complete reference
example.

### What `gpuistl_if_available.hpp` provides

After the include, the following are available inside `#if HAVE_CUDA`
blocks (no further includes needed):

| Symbol | Header it comes from |
| ------ | -------------------- |
| `Opm::gpuistl::GpuBuffer<T>` | `GpuBuffer.hpp` |
| `Opm::gpuistl::GpuView<T>` | `GpuView.hpp` |
| `Opm::gpuistl::GpuVector<T>` | `GpuVector.hpp` |
| `Opm::gpuistl::make_gpu_shared_ptr` / `make_gpu_unique_ptr` | `gpu_smart_pointer.hpp` |
| `Opm::gpuistl::PointerView<T>` / `ValueAsPointer<T>` | `gpu_smart_pointer.hpp` |

---

## 5d. The static/non-static fluid system duality

### Background: why a non-static twin exists

`BlackOilFluidSystem` was historically a **fully static** class: all its
data lives in `static` member variables, and all its methods are
`static`. This design works perfectly on the CPU — there is effectively
a single global instance — but is incompatible with the GPU because
`static` storage cannot be placed in device memory.

To solve this without breaking existing CPU code, a **non-static twin**
`BlackOilFluidSystemNonStatic` was introduced
(`opm-common/opm/material/fluidsystems/BlackOilFluidSystemNonStatic.hpp`).
Both classes share the same implementation via a macro-template
(`BlackOilFluidSystem_macrotemplate.hpp`), so they are always in sync.
The static class exposes `getNonStaticInstance()` which returns a
reference to a `BlackOilFluidSystemNonStatic` object that holds the same
data as the static class but as ordinary member variables — making it
portable to the GPU via the standard `copy_to_gpu` / `make_view` pattern.

The two classes are **interchangeable as template arguments**: any class
that is templated on `FluidSystemT` can be instantiated with either the
static (`BlackOilFluidSystem`) or the non-static
(`BlackOilFluidSystemNonStatic<Scalar, IndexTraits, GpuView>`) variant.

### How dependent classes support both variants: `fluidSystemIsStatic`

`BlackOilFluidState` illustrates the canonical pattern for a class that
needs to call into a fluid system but must work with *both* the static
and the non-static variant:

```cpp
// FluidSystemT may be BlackOilFluidSystem (static) or
// BlackOilFluidSystemNonStatic<..., GpuView> (GPU view).
template <class ValueT, class FluidSystemT, ...>
class BlackOilFluidState
{
    using FluidSystem = FluidSystemT;

    // True when FluidSystem is a zero-sized empty type (i.e. fully static).
    static constexpr bool fluidSystemIsStatic = std::is_empty_v<FluidSystem>;

    // Unified accessor — compiles and works for both variants.
    OPM_HOST_DEVICE const FluidSystem& fluidSystem() const
    {
        if constexpr (fluidSystemIsStatic) {
            static FluidSystem instance;
            return instance;          // zero-cost; no stored pointer
        } else {
            return **fluidSystemPtr_; // double-deref: ConditionalStorage<T*>
        }
    }

    // Only occupies storage when FluidSystem is not a zero-sized empty type.
    ConditionalStorage<!fluidSystemIsStatic, FluidSystem const*> fluidSystemPtr_;
};
```

When instantiated with the **static** fluid system, `fluidSystemIsStatic`
is `true`, `ConditionalStorage` stores nothing, and `fluidSystem()`
returns a `static` local — no overhead, no stored pointer, fully
backwards-compatible.

When instantiated with the **non-static / GPU view** variant,
`fluidSystemIsStatic` is `false`, the constructor stores a pointer to the
provided fluid system, and `fluidSystem()` dereferences it. Because the
GPU-view fluid system is trivially copyable and its data lives in device
memory, the pointer is valid inside a kernel.

### CPU-to-GPU flow for a class with a non-static fluid system

```cpp
// 1. Initialise the static fluid system as usual.
FluidSystem::initFromState(eclState, schedule);

// 2. Obtain the non-static twin — owned by the static fluid system,
//    valid as long as no re-init occurs.
auto& dynamicFs = FluidSystem::getNonStaticInstance();

// 3. Copy the fluid system to the GPU (owning buffer).
auto gpuFsBuffer = Opm::gpuistl::copy_to_gpu(dynamicFs);

// 4. Create a trivially copyable view to pass to kernels.
auto gpuFsView = Opm::gpuistl::make_view(gpuFsBuffer);

// 5. Construct the dependent object using the non-static GPU view type.
//    The fluid-state is now parameterised on the GpuView fluid system.
using GpuFluidSystem = decltype(gpuFsView);
BlackOilFluidState<Scalar, GpuFluidSystem, ...> gpuFluidState(gpuFsView);

// 6. Pass by value to the kernel (gpuFsBuffer must outlive the kernel).
myKernel<<<grid, block>>>(gpuFluidState, /* … */);
OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
```

The complete working example is in
`opm-simulators/tests/gpuistl/test_gpuBlackOilFluidSystem.cu`.

---

## 6. Putting it together: launching a kernel

A complete CPU-side example, mirroring
`tests/gpuistl/test_gpu_linear_two_phase_material.cu`:

```cpp
using Scalar = float;

// 1. Build and finalize the CPU object.
auto cpu = Opm::MyComposite<Traits>({0.0f, 0.5f, 1.0f}, {0.0f, 0.9f, 1.0f});

// 2. Owning GPU copy — no template arg: return type is always MyComposite<Traits, GpuBuffer>.
//    Keep this alive until all kernels are done.
auto gpuOwned = Opm::gpuistl::copy_to_gpu(cpu);

// 3. Cheap, trivially copyable view — no template arg: return type is always MyComposite<Traits, GpuView>.
//    Pass by value to the kernel.
auto gpuView = Opm::gpuistl::make_view(gpuOwned);

// 4. Launch.
Scalar* dResult = nullptr;
OPM_GPU_SAFE_CALL(cudaMalloc(&dResult, sizeof(Scalar)));

evaluateKernel<<<1, 1>>>(gpuView, /*x=*/0.25f, dResult);
OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());

Scalar hResult = 0;
OPM_GPU_SAFE_CALL(cudaMemcpy(&hResult, dResult,
                             sizeof(Scalar), cudaMemcpyDeviceToHost));
OPM_GPU_SAFE_CALL(cudaFree(dResult));

// 5. gpuOwned destroyed here → GPU memory released.
```

The kernel itself simply dereferences the view:

```cpp
using GpuViewComposite = Opm::MyComposite<Traits, Opm::gpuistl::GpuView>;

__global__
void evaluateKernel(GpuViewComposite params, Scalar x, Scalar* out)
{
    *out = params.evaluate(x);
}
```

Always wrap CUDA/HIP runtime calls in **`OPM_GPU_SAFE_CALL`** (from
`opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp`) so that
errors are turned into `OPM_THROW`-style exceptions on the host.

---

## 7. Adding tests

Tests live in **`opm-simulators/tests/gpuistl/`** — that is the only
place where `.cu` files are actually compiled. The pattern is:

1. Create `tests/gpuistl/test_my_composite.cu`. Use Boost.Test as in
   the surrounding files. The header skeleton looks like:

   ```cpp
   #include <config.h>

   #define BOOST_TEST_MODULE TestMyComposite
   #include <boost/test/unit_test.hpp>

   #include <cuda_runtime.h>

   #include <opm/path/to/MyComposite.hpp>
   #include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
   #include <opm/simulators/linalg/gpuistl/GpuView.hpp>
   #include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
   ```

2. Define one or more `__global__` wrapper kernels that take the *view*
   instantiation of your class and write a result into device memory.
3. Inside `BOOST_AUTO_TEST_CASE`, build the CPU object, call
   `copy_to_gpu`, call `make_view`, launch the kernel, copy the result
   back, and `BOOST_CHECK` against the CPU reference computed from the
   exact same inputs. Always compare GPU output to a CPU reference
   computed from the *same* class — that is what guarantees identical
   semantics on both backends.
4. Register the test in
   **`opm-simulators/CMakeLists_files.cmake`** by adding a single
   line next to the existing GPU tests, e.g.:

   ```cmake
   ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_my_composite.cu)
   ```

   The `ADD_CUDA_OR_HIP_FILE` macro takes care of wiring the file into
   the CUDA build *and* into the HIP build (see section 8).

5. Add the test name to the `foreach(test ...)` block inside the
   `opm-simulators_tests_hook` macro in
   **`opm-simulators/CMakeLists.txt`** (around line 512). That loop
   tags every GPU test with the `gpu_cuda` or `gpu_hip` CTest label
   (selected automatically based on `CONVERT_CUDA_TO_HIP`), which is
   what the CI uses to pick up GPU tests. The entry is just the bare
   test name (without the `test_` prefix), inserted in alphabetical
   order, e.g.:

   ```cmake
   foreach(test
                ...
                MiniVector
                my_composite
                preconditioner_factory_gpu
                ...
          )
   ```

   Without this entry the test will still build and run locally, but
   it will not be labelled as a GPU test and therefore will be skipped
   by GPU-targeted CI jobs.

Reference test files to read before writing your own:

* `tests/gpuistl/test_gpu_linear_two_phase_material.cu` — composite
  class with multiple `GpuBuffer`s.
* `tests/gpuistl/test_GpuSparseTable.cu` — composite class living in
  `opm-grid`.
* `tests/gpuistl/test_GpuBuffer.cu` — exercises the buffer/view
  primitives themselves.
* `tests/gpuistl/test_throw_macros_on_gpu.cu` — shows the expected
  behaviour of `OPM_THROW` from device code.

---

## 8. Relation to hipification

When the build is configured for AMD GPUs, the CUDA sources are
translated to HIP at build time. The mechanism lives **entirely inside
`opm-simulators`**:

* `opm-simulators/CMakeLists.txt` looks for `hipify-perl`,
  imports it as the `hipify-perl` CMake target and, when
  `CONVERT_CUDA_TO_HIP` is `ON`, registers the `hipified_headers`
  target that the GPU library/test executables depend on.
* `opm-simulators/CMakeLists_files.cmake` defines the macro
  `ADD_CUDA_OR_HIP_FILE(LIST DIR FILE)`. When CUDA is found the file
  is added as-is; when HIP is the active backend the macro emits an
  `add_custom_command` that runs `bin/hipify_file.sh` to translate
  `${DIR}/gpuistl/${FILE}` into
  `${PROJECT_BINARY_DIR}/${DIR}/gpuistl_hip/${FILE/.cu/.hip}` and adds
  the generated file to the build. Public headers are tracked
  separately and end up as a `hipified_headers` target dependency.
* `opm-simulators/bin/hipify_file.sh` is the actual translation step:
  it runs `hipify-perl` and post-processes the `#include` paths so
  they refer to the `gpuistl_hip/` (rather than `gpuistl/`) versions of
  any sibling files.

Practical implications when porting a class:

* **`opm-common` and `opm-grid` headers are *not* hipified.** They are
  consumed unchanged from the hipified `opm-simulators` translation
  units. This is exactly why `OPM_HOST_DEVICE`, `OPM_THROW`,
  `<cuda_runtime.h>` vs `<hip/hip_runtime.h>` and the `HAVE_CUDA`
  guards are routed through `gpuDecorators.hpp` — those headers must
  compile cleanly under both backends. Do **not** `#include` runtime
  headers directly from a header that lives in `opm-common` or
  `opm-grid`.
* If you add a new GPU source file to `opm-simulators`, register it
  through `ADD_CUDA_OR_HIP_FILE` so that it is picked up automatically
  by the HIP build. Never list `.cu` files unconditionally.
* If you add a new public header that should *not* be hipified (for
  example because it must remain reachable from non-GPU translation
  units that are built without HIP), add it to `PUBLIC_HEADER_FILES`
  directly — see the comment near
  `opm/simulators/linalg/gpuistl/device_management.hpp` in
  `CMakeLists_files.cmake` for the precedent.
* The `OPM_IS_COMPILING_WITH_GPU_COMPILER` macro is `1` under both
  `nvcc` and `hipcc`, so it is the right knob if you ever need to
  branch on "GPU compiler vs host compiler" rather than on
  "CUDA vs HIP".

In short: write the class for CUDA in `opm-simulators` (and pure,
header-only, decorator-guarded code in `opm-common` / `opm-grid`),
register sources through `ADD_CUDA_OR_HIP_FILE`, and the HIP backend
falls out automatically.

---

## Checklist

- [ ] **The existing class file was edited in-place** — no duplicate
      GPU-only class was created. A separate GPU class is only
      acceptable when the original class has a heavily virtual interface
      or is otherwise impossible to template.
- [ ] Class is templated on its storage container using
      `template<class> class Storage = VectorWithDefaultAllocator`.
- [ ] The CPU default instantiation (`Storage = VectorWithDefaultAllocator`)
      is behaviourally identical to the class before the porting work.
- [ ] Only the methods that must run on the device are annotated with
      `OPM_HOST_DEVICE`; host-only methods are left undecorated.
- [ ] Host-only methods that make no sense for a `GpuBuffer` /
      `GpuView` instantiation are guarded (e.g. `static_assert` or
      `if constexpr`) so they fail clearly at compile time if misused.
- [ ] Every member function callable from a kernel is annotated with
      `OPM_HOST_DEVICE`.
- [ ] All exceptions go through `OPM_THROW`.
- [ ] `<opm/common/utility/gpuDecorators.hpp>` is the only GPU-related
      header included from the class definition itself.
- [ ] `copy_to_gpu` and `make_view` live in `namespace Opm::gpuistl`,
      are friended by the class, and are guarded by `#if HAVE_CUDA`
      when the class itself is in `opm-common` or `opm-grid`.
- [ ] GPU headers in `opm-common` / `opm-grid` are pulled in via
      `<opm/common/utility/gpuistl_if_available.hpp>` (unconditionally
      in `opm-common`; guarded by `#if HAVE_OPM_COMMON` in `opm-grid`)
      rather than manually including individual `GpuBuffer` / `GpuView`
      headers inside `#if HAVE_CUDA` / `#if USE_HIP` blocks.
- [ ] The CPU code that calls `copy_to_gpu` keeps the returned owning
      object alive until every kernel using the view has completed.
- [ ] A new `test_*.cu` exists in `opm-simulators/tests/gpuistl/` and
      is registered via `ADD_CUDA_OR_HIP_FILE` in
      `opm-simulators/CMakeLists_files.cmake`.
- [ ] The test name is added to the `foreach(test ...)` block in
      `opm-simulators/CMakeLists.txt` so it gets the `gpu_cuda` /
      `gpu_hip` CTest label.