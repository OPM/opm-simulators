// SPDX-License-Identifier: GPL-3.0-or-later
#include "VulkanSolver.hpp"
#include "PressureHierarchy.hpp"
#include "bicgstab_dilu_spv.hpp"
#include "bicgstab_spv.hpp"
#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vulkan/vulkan.h>

namespace Reslab
{
namespace
{
    void check(VkResult result, const char* operation)
    {
        if (result != VK_SUCCESS)
            throw std::runtime_error(std::string(operation)
                                     + " failed, VkResult=" + std::to_string(result));
    }
    std::string lower(std::string s)
    {
        std::transform(
            s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
        return s;
    }
    using Clock = std::chrono::steady_clock;
} // namespace
struct VulkanSolver::Impl {
    struct Buffer {
        VkBuffer buffer {};
        VkDeviceMemory memory {};
        void* mapped {};
        VkDeviceSize size {};
    };
    VkInstance instance {};
    VkPhysicalDevice physical {};
    VkPhysicalDeviceProperties properties {};
    VkDevice device {};
    VkQueue queue {};
    uint32_t family {};
    VkDescriptorSetLayout setLayout {};
    VkDescriptorPool descriptorPool {};
    VkDescriptorSet descriptor {};
    VkPipelineLayout layout {};
    VkShaderModule shader {};
    std::array<VkPipeline, 30> pipelines {};
    VkCommandPool commandPool {};
    VkCommandBuffer command {};
    VkFence fence {};
    VkQueryPool timestamps {};
    bool profile = false, profiled = false;
    uint32_t timestampBits = 0;
    std::vector<uint32_t> timedOperations;
    std::array<Buffer, 8> buffers {};
    std::vector<uint32_t> rows, columns;
    std::vector<uint32_t> inputRows, inputColumns, permutation, valueOrder;
    std::vector<uint32_t> fullRows, fullColumns;
    std::vector<double> values, originalValues, scale, inverse;
    std::string name;
    uint32_t n = 0;
    bool prepared = false;
    bool dilu = false;
    bool sweeps = false;
    bool cpr = false;
    bool cooperative = false;
    bool fpf = false;
    uint32_t sweepCount = 5, pressureCycles = 1;
    PressureHierarchy pressure;
    std::vector<uint32_t> blockRows, blockCols, scalarToBlock, levelRows, lowerOffsets,
        upperOffsets;
    std::vector<double> blockValues;
    uint32_t bufferCount() const
    {
        return dilu ? 8 : 5;
    }

    ~Impl()
    {
        if (device) {
            vkDeviceWaitIdle(device);
            releaseBuffers();
            if (fence)
                vkDestroyFence(device, fence, nullptr);
            if (timestamps)
                vkDestroyQueryPool(device, timestamps, nullptr);
            if (commandPool)
                vkDestroyCommandPool(device, commandPool, nullptr);
            for (auto pipeline : pipelines)
                if (pipeline)
                    vkDestroyPipeline(device, pipeline, nullptr);
            if (shader)
                vkDestroyShaderModule(device, shader, nullptr);
            if (layout)
                vkDestroyPipelineLayout(device, layout, nullptr);
            if (descriptorPool)
                vkDestroyDescriptorPool(device, descriptorPool, nullptr);
            if (setLayout)
                vkDestroyDescriptorSetLayout(device, setLayout, nullptr);
            vkDestroyDevice(device, nullptr);
        }
        if (instance)
            vkDestroyInstance(instance, nullptr);
    }
    void releaseBuffers()
    {
        for (auto& b : buffers) {
            if (b.mapped)
                vkUnmapMemory(device, b.memory);
            if (b.buffer)
                vkDestroyBuffer(device, b.buffer, nullptr);
            if (b.memory)
                vkFreeMemory(device, b.memory, nullptr);
            b = {};
        }
    }
    void init(const std::string& match)
    {
        VkApplicationInfo app {VK_STRUCTURE_TYPE_APPLICATION_INFO};
        app.pApplicationName = "OPM Flow Vulkan experimental";
        app.apiVersion = VK_API_VERSION_1_1;
        VkInstanceCreateInfo ici {VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
        ici.pApplicationInfo = &app;
        check(vkCreateInstance(&ici, nullptr, &instance), "vkCreateInstance");
        uint32_t count = 0;
        check(vkEnumeratePhysicalDevices(instance, &count, nullptr), "enumerate devices");
        std::vector<VkPhysicalDevice> devices(count);
        check(vkEnumeratePhysicalDevices(instance, &count, devices.data()), "enumerate devices");
        int matches = 0;
        for (auto candidate : devices) {
            VkPhysicalDeviceProperties p;
            vkGetPhysicalDeviceProperties(candidate, &p);
            if (p.deviceType != VK_PHYSICAL_DEVICE_TYPE_CPU
                && lower(p.deviceName).find(lower(match)) != std::string::npos) {
                physical = candidate;
                properties = p;
                ++matches;
            }
        }
        if (matches != 1)
            throw std::runtime_error("Vulkan requires exactly one hardware device matching '"
                                     + match + "'; CPU fallback forbidden");
        name = properties.deviceName;
        VkPhysicalDeviceFeatures features;
        vkGetPhysicalDeviceFeatures(physical, &features);
        if (!features.shaderFloat64)
            throw std::runtime_error("Vulkan device lacks shaderFloat64");
        vkGetPhysicalDeviceQueueFamilyProperties(physical, &count, nullptr);
        std::vector<VkQueueFamilyProperties> families(count);
        vkGetPhysicalDeviceQueueFamilyProperties(physical, &count, families.data());
        family = count;
        for (uint32_t i = 0; i < count; ++i)
            if (families[i].queueFlags & VK_QUEUE_COMPUTE_BIT) {
                family = i;
                break;
            }
        if (family == count)
            throw std::runtime_error("No Vulkan compute queue");
        timestampBits = families[family].timestampValidBits;
        profile = std::getenv("RESLAB_VULKAN_PROFILE") != nullptr;
        if (profile && !timestampBits)
            throw std::runtime_error("Vulkan queue has no timestamps");
        float priority = 1;
        VkDeviceQueueCreateInfo qci {VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
        qci.queueFamilyIndex = family;
        qci.queueCount = 1;
        qci.pQueuePriorities = &priority;
        VkPhysicalDeviceFeatures enabled {};
        enabled.shaderFloat64 = VK_TRUE;
        VkDeviceCreateInfo dci {VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
        dci.queueCreateInfoCount = 1;
        dci.pQueueCreateInfos = &qci;
        dci.pEnabledFeatures = &enabled;
        check(vkCreateDevice(physical, &dci, nullptr, &device), "vkCreateDevice");
        vkGetDeviceQueue(device, family, 0, &queue);
        std::array<VkDescriptorSetLayoutBinding, 8> bindings {};
        for (uint32_t i = 0; i < bufferCount(); ++i)
            bindings[i]
                = {i, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};
        VkDescriptorSetLayoutCreateInfo sl {VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
        sl.bindingCount = bufferCount();
        sl.pBindings = bindings.data();
        check(vkCreateDescriptorSetLayout(device, &sl, nullptr, &setLayout), "descriptor layout");
        VkDescriptorPoolSize poolSize {VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, bufferCount()};
        VkDescriptorPoolCreateInfo dp {VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
        dp.maxSets = 1;
        dp.poolSizeCount = 1;
        dp.pPoolSizes = &poolSize;
        check(vkCreateDescriptorPool(device, &dp, nullptr, &descriptorPool), "descriptor pool");
        VkDescriptorSetAllocateInfo da {VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
        da.descriptorPool = descriptorPool;
        da.descriptorSetCount = 1;
        da.pSetLayouts = &setLayout;
        check(vkAllocateDescriptorSets(device, &da, &descriptor), "allocate descriptors");
        VkPushConstantRange range {VK_SHADER_STAGE_COMPUTE_BIT, 0, dilu ? 16u : 8u};
        VkPipelineLayoutCreateInfo pl {VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        pl.setLayoutCount = 1;
        pl.pSetLayouts = &setLayout;
        pl.pushConstantRangeCount = 1;
        pl.pPushConstantRanges = &range;
        check(vkCreatePipelineLayout(device, &pl, nullptr, &layout), "pipeline layout");
        VkShaderModuleCreateInfo sm {VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
        sm.codeSize = dilu ? sizeof(bicgstab_dilu_spv) : sizeof(bicgstab_spv);
        sm.pCode = dilu ? bicgstab_dilu_spv : bicgstab_spv;
        check(vkCreateShaderModule(device, &sm, nullptr, &shader), "shader module");
        VkComputePipelineCreateInfo pc {VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
        pc.layout = layout;
        pc.stage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        pc.stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
        pc.stage.module = shader;
        pc.stage.pName = "main";
        cooperative = std::getenv("RESLAB_VULKAN_COOPERATIVE") != nullptr;
        fpf = std::getenv("RESLAB_VULKAN_CPR_FPF") != nullptr;
        if (const char* s = std::getenv("RESLAB_VULKAN_SWEEP_COUNT")) {
            std::string value(s);
            if (value != "1" && value != "3" && value != "5" && value != "7" && value != "9")
                throw std::invalid_argument("Sweep count must be 1, 3, 5, 7 or 9");
            sweepCount = std::stoul(value);
        }
        if (const char* s = std::getenv("RESLAB_VULKAN_PRESSURE_CYCLES")) {
            std::string value(s);
            if (value != "1" && value != "2" && value != "3")
                throw std::invalid_argument("Pressure cycles must be 1, 2 or 3");
            pressureCycles = std::stoul(value);
        }
        if (fpf && !cpr)
            throw std::invalid_argument("FPF requires CPR");
        if (cpr)
            std::clog << "[Vulkan] CPR stages=" << (fpf ? "fine-pressure-fine" : "pressure-fine")
                      << '\n';
        VkSpecializationMapEntry entries[3] = {{0, 0, sizeof(uint32_t)},
                                               {1, sizeof(uint32_t), sizeof(uint32_t)},
                                               {2, 2 * sizeof(uint32_t), sizeof(uint32_t)}};
        uint32_t constants[3] = {
            0, std::getenv("RESLAB_VULKAN_PRECONDITIONER_FP32") ? 1u : 0u, cooperative ? 1u : 0u};
        if (cooperative && constants[1])
            throw std::invalid_argument("Cooperative sweeps require FP64");
        std::clog << "[Vulkan] DILU_sweeps_arithmetic=" << (constants[1] ? "FP32" : "FP64")
                  << "; outer_solver=FP64; true_residual=FP64\n";
        VkSpecializationInfo specialization {};
        specialization.mapEntryCount = 3;
        specialization.pMapEntries = entries;
        specialization.dataSize = sizeof(constants);
        specialization.pData = constants;
        for (uint32_t op = 0; op < (dilu ? 30u : 13u); ++op) {
            const bool used = (op < 13 && !(dilu && (op == 1 || op == 6)))
                || (op >= 13 && op <= 16 && dilu && !sweeps) || (op >= 17 && op <= 20 && sweeps)
                || (op >= 21 && op <= 28 && cpr) || (op == 29 && fpf);
            if (!used)
                continue;
            constants[0] = op;
            pc.stage.pSpecializationInfo = &specialization;
            check(vkCreateComputePipelines(device, VK_NULL_HANDLE, 1, &pc, nullptr, &pipelines[op]),
                  "specialized compute pipeline");
        }
        VkCommandPoolCreateInfo cp {VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
        cp.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
        cp.queueFamilyIndex = family;
        check(vkCreateCommandPool(device, &cp, nullptr, &commandPool), "command pool");
        VkCommandBufferAllocateInfo ca {VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
        ca.commandPool = commandPool;
        ca.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
        ca.commandBufferCount = 1;
        check(vkAllocateCommandBuffers(device, &ca, &command), "command buffer");
        VkFenceCreateInfo fc {VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
        check(vkCreateFence(device, &fc, nullptr, &fence), "fence");
        if (profile) {
            VkQueryPoolCreateInfo qi {VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO};
            qi.queryType = VK_QUERY_TYPE_TIMESTAMP;
            qi.queryCount = 32768;
            check(vkCreateQueryPool(device, &qi, nullptr, &timestamps), "timestamp pool");
        }
    }
    void allocate(Buffer& b, size_t bytes)
    {
        b.size = std::max<size_t>(bytes, 8);
        if (b.size > properties.limits.maxStorageBufferRange)
            throw std::runtime_error("Vulkan storage buffer size limit exceeded");
        VkBufferCreateInfo bi {VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
        bi.size = b.size;
        bi.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT;
        bi.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
        check(vkCreateBuffer(device, &bi, nullptr, &b.buffer), "buffer");
        VkMemoryRequirements req;
        vkGetBufferMemoryRequirements(device, b.buffer, &req);
        VkPhysicalDeviceMemoryProperties memory;
        vkGetPhysicalDeviceMemoryProperties(physical, &memory);
        const auto flags
            = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;
        uint32_t index = memory.memoryTypeCount;
        for (uint32_t i = 0; i < memory.memoryTypeCount; ++i) {
            if ((req.memoryTypeBits & (1u << i))
                && (memory.memoryTypes[i].propertyFlags & flags) == flags) {
                index = i;
                if (memory.memoryTypes[i].propertyFlags & VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT)
                    break;
            }
        }
        if (index == memory.memoryTypeCount)
            throw std::runtime_error("No host coherent Vulkan storage memory");
        VkMemoryAllocateInfo ma {VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
        ma.allocationSize = req.size;
        ma.memoryTypeIndex = index;
        check(vkAllocateMemory(device, &ma, nullptr, &b.memory), "allocate memory");
        check(vkBindBufferMemory(device, b.buffer, b.memory, 0), "bind memory");
        check(vkMapMemory(device, b.memory, 0, b.size, 0, &b.mapped), "map memory");
    }
    void barrier(VkPipelineStageFlags src,
                 VkPipelineStageFlags dst,
                 VkAccessFlags srcAccess,
                 VkAccessFlags dstAccess)
    {
        VkMemoryBarrier b {VK_STRUCTURE_TYPE_MEMORY_BARRIER};
        b.srcAccessMask = srcAccess;
        b.dstAccessMask = dstAccess;
        vkCmdPipelineBarrier(command, src, dst, 0, 1, &b, 0, nullptr, 0, nullptr);
    }
    void record()
    {
        check(vkResetCommandBuffer(command, 0), "reset command buffer");
        VkCommandBufferBeginInfo begin {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
        check(vkBeginCommandBuffer(command, &begin), "begin command buffer");
        timedOperations.clear();
        if (profile)
            vkCmdResetQueryPool(command, timestamps, 0, 32768);
        barrier(VK_PIPELINE_STAGE_HOST_BIT | VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                VK_ACCESS_HOST_WRITE_BIT | VK_ACCESS_SHADER_WRITE_BIT,
                VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT);
        vkCmdBindDescriptorSets(
            command, VK_PIPELINE_BIND_POINT_COMPUTE, layout, 0, 1, &descriptor, 0, nullptr);
        auto dispatch = [&](uint32_t op, uint32_t offset, uint32_t count, uint32_t groups) {
            if (groups > properties.limits.maxComputeWorkGroupCount[0])
                throw std::runtime_error("Vulkan dispatch size limit exceeded");
            if (profile) {
                if (timedOperations.size() >= 16384)
                    throw std::runtime_error("Too many profiled dispatches");
                vkCmdWriteTimestamp(command,
                                    VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                                    timestamps,
                                    timedOperations.size() * 2);
            }
            vkCmdBindPipeline(command, VK_PIPELINE_BIND_POINT_COMPUTE, pipelines[op]);
            uint32_t push[4] = {n, op, offset, count};
            vkCmdPushConstants(
                command, layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, dilu ? 16 : 8, push);
            vkCmdDispatch(command, groups, 1, 1);
            if (profile) {
                vkCmdWriteTimestamp(command,
                                    VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                                    timestamps,
                                    timedOperations.size() * 2 + 1);
                timedOperations.push_back(op);
            }
            barrier(VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                    VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                    VK_ACCESS_SHADER_WRITE_BIT,
                    VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT);
        };
        for (unsigned iteration = 0; iteration < 8; ++iteration)
            for (uint32_t op = 0; op < 13; ++op) {
                if (dilu && (op == 1 || op == 6)) {
                    if (cpr) {
                        if (fpf) {
                            uint32_t cellsPerGroup = cooperative ? 8 : 64;
                            for (uint32_t k = 0; k < sweepCount; ++k)
                                dispatch(op == 1 ? 17 : 19,
                                         k,
                                         n / 3,
                                         (n / 3 + cellsPerGroup - 1) / cellsPerGroup);
                            for (uint32_t k = 0; k + 1 < sweepCount; ++k)
                                dispatch(op == 1 ? 18 : 20,
                                         k,
                                         n / 3,
                                         (n / 3 + cellsPerGroup - 1) / cellsPerGroup);
                            dispatch(29, 0, op, (n + 63) / 64);
                        }
                        dispatch(21, 0, op + (fpf ? 32 : 0), (n / 3 + 63) / 64);
                        for (uint32_t cycle = 0; cycle < pressureCycles; ++cycle) {
                            for (size_t l = 0; l + 1 < pressure.levels.size(); ++l) {
                                auto& level = pressure.levels[l];
                                uint32_t groups = (level.size() + 63) / 64;
                                dispatch(22, level.meta, op, groups);
                                dispatch(27, level.meta, op, groups);
                                dispatch(
                                    23, level.meta, op, (pressure.levels[l + 1].size() + 63) / 64);
                            }
                            auto& coarse = pressure.levels.back();
                            dispatch(25, coarse.meta, op, (coarse.size() + 63) / 64);
                            for (size_t l = pressure.levels.size() - 1; l-- > 0;) {
                                auto& level = pressure.levels[l];
                                uint32_t groups = (level.size() + 63) / 64;
                                dispatch(24, level.meta, op, groups);
                                dispatch(22, level.meta, op, groups);
                                dispatch(27, level.meta, op, groups);
                            }
                        }
                        dispatch(26, 0, op + (fpf ? 32 : 0), (n + 63) / 64);
                    }
                    if (sweeps) {
                        uint32_t cellsPerGroup = cooperative ? 8 : 64;
                        for (uint32_t k = 0; k < sweepCount; ++k)
                            dispatch(op == 1 ? 17 : 19,
                                     k,
                                     cpr ? 0 : n / 3,
                                     (n / 3 + cellsPerGroup - 1) / cellsPerGroup);
                        for (uint32_t k = 0; k + 1 < sweepCount; ++k)
                            dispatch(op == 1 ? 18 : 20,
                                     k,
                                     cpr ? 0 : n / 3,
                                     (n / 3 + cellsPerGroup - 1) / cellsPerGroup);
                        if (cpr)
                            dispatch(28, 0, op + (fpf ? 32 : 0), (n / 3 + 63) / 64);
                        continue;
                    }
                    for (size_t k = 0; k + 1 < lowerOffsets.size(); ++k) {
                        uint32_t count = lowerOffsets[k + 1] - lowerOffsets[k];
                        dispatch(op == 1 ? 13 : 15, lowerOffsets[k], count, (count + 63) / 64);
                    }
                    for (size_t k = 0; k + 1 < upperOffsets.size(); ++k) {
                        uint32_t count = upperOffsets[k + 1] - upperOffsets[k];
                        dispatch(
                            op == 1 ? 14 : 16, n / 3 + upperOffsets[k], count, (count + 63) / 64);
                    }
                } else
                    dispatch(op, 0, 0, (op == 4 || op == 9 || op == 12) ? 1 : (n + 63) / 64);
            }
        barrier(VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                VK_PIPELINE_STAGE_HOST_BIT,
                VK_ACCESS_SHADER_WRITE_BIT,
                VK_ACCESS_HOST_READ_BIT);
        check(vkEndCommandBuffer(command), "end command buffer");
    }
    void resize()
    {
        check(vkDeviceWaitIdle(device), "wait before resizing");
        // Small changes in sparsity must not recreate every Vulkan allocation.
        // Keep buffers with enough capacity; grow only the buffers that need it.
        auto ensure = [&](Buffer& b, size_t bytes) {
            bytes = std::max<size_t>(bytes, 8);
            if (b.buffer && b.size >= bytes)
                return;
            if (bytes > properties.limits.maxStorageBufferRange)
                throw std::runtime_error("Vulkan storage buffer size limit exceeded");
            if (b.mapped)
                vkUnmapMemory(device, b.memory);
            if (b.buffer)
                vkDestroyBuffer(device, b.buffer, nullptr);
            if (b.memory)
                vkFreeMemory(device, b.memory, nullptr);
            b = {};
            size_t capacity = std::min<size_t>(properties.limits.maxStorageBufferRange,
                                               (bytes + bytes / 4 + 255) & ~size_t(255));
            allocate(b, capacity);
        };
        ensure(buffers[0], rows.size() * 4);
        ensure(buffers[1], columns.size() * 4);
        ensure(buffers[2], values.size() * 8);
        ensure(buffers[3], size_t(n) * (fpf ? 17 : 16) * 8);
        ensure(buffers[4], (16 + size_t(3) * ((n + 63) / 64)) * 8);
        if (dilu)
            ensure(buffers[5], levelRows.size() * 4);
        if (dilu) {
            ensure(buffers[6], pressure.indices.size() * 4);
            ensure(buffers[7], pressure.data.size() * 8);
        }
        std::array<VkDescriptorBufferInfo, 8> info {};
        std::array<VkWriteDescriptorSet, 8> writes {};
        for (uint32_t i = 0; i < bufferCount(); ++i) {
            info[i] = {buffers[i].buffer, 0, buffers[i].size};
            writes[i].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            writes[i].dstSet = descriptor;
            writes[i].dstBinding = i;
            writes[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
            writes[i].descriptorCount = 1;
            writes[i].pBufferInfo = &info[i];
        }
        vkUpdateDescriptorSets(device, bufferCount(), writes.data(), 0, nullptr);
        record();
    }
    void setupBlockPattern()
    {
        const uint32_t nb = n / 3;
        blockRows.assign(1, 0);
        blockCols.clear();
        scalarToBlock.resize(values.size());
        std::vector<uint32_t> reverse(nb);
        for (uint32_t i = 0; i < nb; ++i)
            reverse[permutation[i]] = i;
        std::vector<uint32_t> neighbors;
        for (uint32_t cell = 0; cell < nb; ++cell) {
            neighbors.clear();
            uint32_t original = permutation[cell];
            for (uint32_t r = 3 * original; r < 3 * original + 3; ++r)
                for (uint32_t j = fullRows[r]; j < fullRows[r + 1]; ++j)
                    neighbors.push_back(reverse[fullColumns[j] / 3]);
            std::sort(neighbors.begin(), neighbors.end());
            neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
            blockCols.insert(blockCols.end(), neighbors.begin(), neighbors.end());
            blockRows.push_back(blockCols.size());
        }
        mapScalarToBlocks();
        blockValues.resize(blockCols.size() * 9);
        std::vector<uint32_t> lower(nb, 0), upper(nb, 0);
        for (uint32_t i = 0; i < nb; ++i)
            for (uint32_t k = blockRows[i]; k < blockRows[i + 1]; ++k)
                if (blockCols[k] < i)
                    lower[i] = std::max(lower[i], lower[blockCols[k]] + 1);
        for (uint32_t i = nb; i-- > 0;)
            for (uint32_t k = blockRows[i]; k < blockRows[i + 1]; ++k)
                if (blockCols[k] > i)
                    upper[i] = std::max(upper[i], upper[blockCols[k]] + 1);
        levelRows.resize(2 * nb);
        auto flatten = [&](const auto& levels, auto& offsets, uint32_t base) {
            offsets.assign(*std::max_element(levels.begin(), levels.end()) + 2, 0);
            for (auto level : levels)
                ++offsets[level + 1];
            std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
            auto next = offsets;
            for (uint32_t i = 0; i < nb; ++i)
                levelRows[base + next[levels[i]]++] = i;
        };
        flatten(lower, lowerOffsets, 0);
        flatten(upper, upperOffsets, nb);
        std::clog << "[Vulkan] DILU lower_levels=" << lowerOffsets.size() - 1
                  << " upper_levels=" << upperOffsets.size() - 1 << '\n';
    }
    void mapScalarToBlocks()
    {
        scalarToBlock.resize(values.size());
        for (uint32_t cell = 0; cell < n / 3; ++cell)
            for (uint32_t r = 0; r < 3; ++r)
                for (uint32_t j = rows[3 * cell + r]; j < rows[3 * cell + r + 1]; ++j) {
                    auto begin = blockCols.begin() + blockRows[cell],
                         end = blockCols.begin() + blockRows[cell + 1];
                    auto entry = std::lower_bound(begin, end, columns[j] / 3);
                    if (entry == end || *entry != columns[j] / 3)
                        throw std::logic_error("Scalar entry missing from symbolic block graph");
                    scalarToBlock[j] = 9 * (entry - blockCols.begin()) + 3 * r + columns[j] % 3;
                }
    }
    void reorder(const std::vector<uint32_t>& sourceRows, const std::vector<uint32_t>& sourceCols)
    {
        const uint32_t nb = n / 3;
        if (!sweeps || permutation.size() != nb) {
            std::vector<std::vector<uint32_t>> neighbors(nb);
            for (uint32_t i = 0; i < n; ++i)
                for (uint32_t k = sourceRows[i]; k < sourceRows[i + 1]; ++k) {
                    uint32_t a = i / 3, b = sourceCols[k] / 3;
                    if (a != b) {
                        neighbors[a].push_back(b);
                        neighbors[b].push_back(a);
                    }
                }
            std::vector<uint32_t> color(nb, UINT32_MAX);
            uint32_t colors = 0;
            for (uint32_t i = 0; i < nb; ++i) {
                auto& adjacent = neighbors[i];
                std::sort(adjacent.begin(), adjacent.end());
                adjacent.erase(std::unique(adjacent.begin(), adjacent.end()), adjacent.end());
                std::vector<bool> used(adjacent.size() + 1, false);
                for (auto j : adjacent)
                    if (color[j] < used.size())
                        used[color[j]] = true;
                color[i] = std::find(used.begin(), used.end(), false) - used.begin();
                colors = std::max(colors, color[i] + 1);
            }
            permutation.resize(nb);
            std::iota(permutation.begin(), permutation.end(), 0);
            std::stable_sort(permutation.begin(), permutation.end(), [&](auto a, auto b) {
                return color[a] < color[b];
            });
            std::clog << "[Vulkan] DILU graph_colors=" << colors << '\n';
        }
        std::vector<uint32_t> reverse(nb);
        for (uint32_t i = 0; i < nb; ++i)
            reverse[permutation[i]] = i;
        rows.assign(1, 0);
        columns.clear();
        valueOrder.clear();
        rows.reserve(n + 1);
        columns.reserve(sourceCols.size());
        valueOrder.reserve(sourceCols.size());
        std::vector<std::pair<uint32_t, uint32_t>> entries;
        for (uint32_t cell = 0; cell < nb; ++cell)
            for (uint32_t r = 0; r < 3; ++r) {
                const uint32_t original = 3 * permutation[cell] + r;
                entries.clear();
                for (uint32_t k = sourceRows[original]; k < sourceRows[original + 1]; ++k)
                    entries.emplace_back(3 * reverse[sourceCols[k] / 3] + sourceCols[k] % 3, k);
                std::sort(entries.begin(), entries.end());
                for (auto [column, k] : entries) {
                    columns.push_back(column);
                    valueOrder.push_back(k);
                }
                rows.push_back(columns.size());
            }
    }
    void invertBlocks()
    {
        if (dilu) {
            std::fill(blockValues.begin(), blockValues.end(), 0);
            for (size_t k = 0; k < values.size(); ++k)
                blockValues[scalarToBlock[k]] = values[k];
        }
        inverse.assign(size_t(n) * 3, 0);
        for (uint32_t cell = 0; cell < n / 3; ++cell) {
            double a[3][6] {};
            for (uint32_t i = 0; i < 3; ++i) {
                for (uint32_t j = rows[3 * cell + i]; j < rows[3 * cell + i + 1]; ++j)
                    if (columns[j] / 3 == cell)
                        a[i][columns[j] % 3] = values[j];
                a[i][3 + i] = 1;
            }
            if (dilu)
                for (uint32_t k = blockRows[cell]; k < blockRows[cell + 1] && blockCols[k] < cell;
                     ++k) {
                    const uint32_t j = blockCols[k];
                    auto begin = blockCols.begin() + blockRows[j],
                         end = blockCols.begin() + blockRows[j + 1];
                    auto transpose = std::lower_bound(begin, end, cell);
                    if (transpose == end || *transpose != cell)
                        continue;
                    const double* aij = blockValues.data() + 9 * k;
                    const double* aji = blockValues.data() + 9 * (transpose - blockCols.begin());
                    const double* dj = inverse.data() + 9 * j;
                    double tmp[9] {};
                    for (unsigned r = 0; r < 3; ++r)
                        for (unsigned c = 0; c < 3; ++c)
                            for (unsigned t = 0; t < 3; ++t)
                                tmp[3 * r + c] += aij[3 * r + t] * dj[3 * t + c];
                    for (unsigned r = 0; r < 3; ++r)
                        for (unsigned c = 0; c < 3; ++c)
                            for (unsigned t = 0; t < 3; ++t)
                                a[r][c] -= tmp[3 * r + t] * aji[3 * t + c];
                }
            for (unsigned col = 0; col < 3; ++col) {
                unsigned pivot = col;
                for (unsigned row = col + 1; row < 3; ++row)
                    if (std::abs(a[row][col]) > std::abs(a[pivot][col]))
                        pivot = row;
                if (a[pivot][col] == 0 || !std::isfinite(a[pivot][col]))
                    throw std::runtime_error("Singular/nonfinite Vulkan Jacobi diagonal block");
                for (unsigned j = 0; j < 6; ++j)
                    std::swap(a[pivot][j], a[col][j]);
                double divisor = a[col][col];
                for (unsigned j = 0; j < 6; ++j)
                    a[col][j] /= divisor;
                for (unsigned row = 0; row < 3; ++row)
                    if (row != col) {
                        double factor = a[row][col];
                        for (unsigned j = 0; j < 6; ++j)
                            a[row][j] -= factor * a[col][j];
                    }
            }
            for (unsigned i = 0; i < 3; ++i)
                for (unsigned j = 0; j < 3; ++j) {
                    if (!std::isfinite(a[i][j + 3]))
                        throw std::runtime_error("Nonfinite inverse diagonal block");
                    inverse[9 * cell + 3 * i + j] = a[i][j + 3];
                }
        }
    }
};
VulkanSolver::VulkanSolver(const std::string& match, bool useDilu, bool parallelSweeps, bool useCpr)
    : impl_(std::make_unique<Impl>())
{
    if (parallelSweeps && !useDilu)
        throw std::invalid_argument("Parallel triangular sweeps require DILU");
    if (useCpr && !parallelSweeps)
        throw std::invalid_argument("CPR requires parallel DILU sweeps");
    impl_->dilu = useDilu;
    impl_->sweeps = parallelSweeps;
    impl_->cpr = useCpr;
    impl_->init(match);
}
VulkanSolver::~VulkanSolver() = default;
const std::string&
VulkanSolver::deviceName() const
{
    return impl_->name;
}
void
VulkanSolver::prepare(std::vector<uint32_t> rows,
                      std::vector<uint32_t> columns,
                      std::vector<double> values,
                      const std::vector<double>& pressureWeights)
{
    const auto setupStart = Clock::now();
    auto& p = *impl_;
    const bool wasPrepared = p.prepared;
    p.prepared = false;
    if (rows.size() < 4 || (rows.size() - 1) % 3 || rows.size() - 1 > UINT32_MAX / 17
        || values.size() > UINT32_MAX || columns.size() != values.size() || rows.front() != 0
        || rows.back() != values.size())
        throw std::invalid_argument("Invalid Vulkan CSR size (requires 3 equations per cell)");
    if (!std::is_sorted(rows.begin(), rows.end()))
        throw std::invalid_argument("Invalid CSR offsets");
    for (size_t i = 0; i + 1 < rows.size(); ++i) {
        if (rows[i] > rows[i + 1])
            throw std::invalid_argument("Invalid CSR offsets");
        for (uint32_t j = rows[i]; j < rows[i + 1]; ++j)
            if (columns[j] >= rows.size() - 1 || !std::isfinite(values[j])
                || (j > rows[i] && columns[j - 1] >= columns[j]))
                throw std::invalid_argument("CSR requires finite values, sorted unique columns");
    }
    if (!pressureWeights.empty()
        && (pressureWeights.size() != rows.size() - 1
            || !std::all_of(pressureWeights.begin(), pressureWeights.end(), [](double w) {
                   return std::isfinite(w);
               })))
        throw std::invalid_argument("Invalid pressure weights");
    // OPM stores full 3x3 blocks, including structural zeros. Do not send those
    // zeros through every sparse product and triangular sweep.
    const auto validated = Clock::now();
    const bool fullChanged = !wasPrepared || rows != p.fullRows || columns != p.fullColumns;
    if (fullChanged) {
        p.fullRows = rows;
        p.fullColumns = columns;
    }
    size_t nonzero = 0;
    uint32_t begin = 0;
    bool sameSize = p.inputRows.size() == rows.size();
    for (size_t r = 0; r + 1 < rows.size(); ++r) {
        uint32_t end = rows[r + 1];
        rows[r] = nonzero;
        uint32_t old = sameSize ? p.inputRows[r] : 0, oldEnd = sameSize ? p.inputRows[r + 1] : 0;
        for (uint32_t j = begin; j < end; ++j) {
            while (old < oldEnd && p.inputColumns[old] < columns[j])
                ++old;
            // Retain previously used slots when they become zero; coefficients
            // becoming nonzero still expand the pattern and invalidate caches.
            if (values[j] != 0 || (old < oldEnd && p.inputColumns[old] == columns[j])) {
                columns[nonzero] = columns[j];
                values[nonzero++] = values[j];
            }
        }
        begin = end;
    }
    rows.back() = nonzero;
    columns.resize(nonzero);
    values.resize(nonzero);
    const bool changed = !wasPrepared || rows != p.inputRows || columns != p.inputColumns;
    p.n = rows.size() - 1;
    if (changed) {
        p.inputRows = rows;
        p.inputColumns = columns;
        if (p.dilu)
            p.reorder(rows, columns);
        else {
            p.rows = std::move(rows);
            p.columns = std::move(columns);
        }
    }
    if (p.dilu) {
        p.values.resize(values.size());
        for (size_t k = 0; k < values.size(); ++k)
            p.values[k] = values[p.valueOrder[k]];
    } else
        p.values = std::move(values);
    p.originalValues = p.values;
    const auto reordered = Clock::now();
    p.scale.resize(p.n);
    for (uint32_t i = 0; i < p.n; ++i) {
        double maximum = 0;
        for (uint32_t j = p.rows[i]; j < p.rows[i + 1]; ++j)
            maximum = std::max(maximum, std::abs(p.values[j]));
        if (maximum == 0)
            throw std::invalid_argument("Zero matrix row in Vulkan solver");
        p.scale[i] = 1 / maximum;
        for (uint32_t j = p.rows[i]; j < p.rows[i + 1]; ++j)
            p.values[j] *= p.scale[i];
    }
    const auto equilibrated = Clock::now();
    if (p.dilu) {
        if (fullChanged || (changed && !p.sweeps))
            p.setupBlockPattern();
        else if (changed)
            p.mapScalarToBlocks();
    }
    p.invertBlocks();
    const auto factored = Clock::now();
    if (p.cpr) {
        std::vector<double> weights(pressureWeights.size());
        for (size_t i = 0; i < weights.size(); ++i)
            weights[i] = pressureWeights[3 * p.permutation[i / 3] + i % 3] / p.scale[i];
        p.pressure.build(p.rows, p.columns, p.values, weights, p.blockRows, p.blockCols);
        if (changed)
            std::clog << "[Vulkan] pressure_levels=" << p.pressure.levels.size()
                      << " scalar_nonzeros=" << p.values.size()
                      << " weights=" << (weights.empty() ? "quasi-IMPES" : "true-IMPES") << '\n';
    }
    const auto pressureBuilt = Clock::now();
    if (changed || (p.cpr && p.pressure.rebuilt))
        p.resize();
    const auto buffersReady = Clock::now();
    std::memcpy(p.buffers[0].mapped, p.rows.data(), p.rows.size() * 4);
    std::memcpy(p.buffers[1].mapped, p.columns.data(), p.columns.size() * 4);
    std::memcpy(p.buffers[2].mapped, p.values.data(), p.values.size() * 8);
    if (p.dilu)
        std::memcpy(p.buffers[5].mapped, p.levelRows.data(), p.levelRows.size() * 4);
    if (p.cpr) {
        std::memcpy(p.buffers[6].mapped, p.pressure.indices.data(), p.pressure.indices.size() * 4);
        std::memcpy(p.buffers[7].mapped, p.pressure.data.data(), p.pressure.data.size() * 8);
    }
    auto* work = static_cast<double*>(p.buffers[3].mapped);
    for (uint32_t i = 0; i < p.n; ++i)
        work[15 * p.n + i] = 1 / p.scale[i];
    p.prepared = true;
    if (std::getenv("RESLAB_VULKAN_PROFILE_SETUP")) {
        auto ms = [](auto a, auto b) {
            return std::chrono::duration<double, std::milli>(b - a).count();
        };
        std::clog << "[Vulkan setup] validate_ms=" << ms(setupStart, validated)
                  << " reorder_ms=" << ms(validated, reordered)
                  << " scale_ms=" << ms(reordered, equilibrated)
                  << " factor_ms=" << ms(equilibrated, factored)
                  << " pressure_ms=" << ms(factored, pressureBuilt)
                  << " buffers_ms=" << ms(pressureBuilt, buffersReady)
                  << " upload_ms=" << ms(buffersReady, Clock::now()) << '\n';
    }
}
SolveResult
VulkanSolver::solve(const std::vector<double>& rhs,
                    std::vector<double>& x,
                    double tolerance,
                    int maxIterations)
{
    auto& p = *impl_;
    if (!p.prepared || !p.n || rhs.size() != p.n || !(tolerance > 0 && tolerance < 1)
        || maxIterations < 1)
        throw std::invalid_argument("Invalid Vulkan solve arguments");
    const auto start = Clock::now();
    double norm2 = 0;
    for (double b : rhs) {
        if (!std::isfinite(b))
            throw std::invalid_argument("Nonfinite RHS");
        norm2 += b * b;
    }
    if (!std::isfinite(norm2)
        || (norm2 == 0 && std::any_of(rhs.begin(), rhs.end(), [](double b) { return b != 0; })))
        throw std::invalid_argument("RHS norm overflows/underflows FP64");
    auto* w = static_cast<double*>(p.buffers[3].mapped);
    auto* z = static_cast<double*>(p.buffers[4].mapped);
    VkSubmitInfo submit {VK_STRUCTURE_TYPE_SUBMIT_INFO};
    submit.commandBufferCount = 1;
    submit.pCommandBuffers = &p.command;
    SolveResult result;
    x.assign(p.n, 0);
    std::vector<double> solveRhs = rhs;
    if (p.dilu)
        for (uint32_t i = 0; i < p.n / 3; ++i)
            for (uint32_t c = 0; c < 3; ++c)
                solveRhs[3 * i + c] = rhs[3 * p.permutation[i] + c];
    std::vector<double> residual = solveRhs;
    // Equilibrate rows without changing OPM's matrix. Verify the original system;
    // restart on its true residual if the recurrent residual drifts. All passes
    // share the caller's iteration budget, and no CPU linear solver is used.
    for (int restart = 0; restart < 5 && result.iterations < maxIterations; ++restart) {
        std::fill_n(w, size_t(p.n) * 9, 0);
        std::copy(p.inverse.begin(), p.inverse.end(), w + 9 * p.n);
        double scaledNorm2 = 0;
        for (uint32_t i = 0; i < p.n; ++i) {
            w[p.n + i] = residual[i] * p.scale[i];
            w[2 * p.n + i] = w[p.n + i];
            scaledNorm2 += w[p.n + i] * w[p.n + i];
        }
        if (!std::isfinite(scaledNorm2))
            throw std::runtime_error("Scaled RHS norm overflows FP64");
        std::fill_n(z, p.buffers[4].size / 8, 0);
        const int remaining = maxIterations - result.iterations;
        // Krylov arithmetic uses the equilibrated system, but stopping uses
        // the original residual norm. Corrections share the original absolute
        // target instead of demanding another full relative reduction.
        z[0] = scaledNorm2;
        z[1] = 1;
        z[2] = 1;
        z[4] = tolerance * tolerance * norm2;
        z[5] = scaledNorm2 == 0 ? 1 : 0;
        z[7] = remaining;
        int batches = 0;
        while (z[5] == 0 && batches < (remaining + 7) / 8) {
            check(vkResetFences(p.device, 1, &p.fence), "reset fence");
            check(vkQueueSubmit(p.queue, 1, &submit, p.fence), "submit solve");
            check(vkWaitForFences(p.device, 1, &p.fence, VK_TRUE, 30000000000ULL), "wait solve");
            if (p.profile && !p.profiled) {
                std::vector<uint64_t> ticks(p.timedOperations.size() * 2);
                check(vkGetQueryPoolResults(p.device,
                                            p.timestamps,
                                            0,
                                            ticks.size(),
                                            ticks.size() * 8,
                                            ticks.data(),
                                            8,
                                            VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT),
                      "read kernel timestamps");
                std::array<double, 30> milliseconds {};
                uint64_t mask
                    = p.timestampBits == 64 ? UINT64_MAX : (uint64_t(1) << p.timestampBits) - 1;
                for (size_t k = 0; k < p.timedOperations.size(); ++k)
                    milliseconds[p.timedOperations[k]] += ((ticks[2 * k + 1] - ticks[2 * k]) & mask)
                        * p.properties.limits.timestampPeriod * 1e-6;
                for (size_t op = 0; op < milliseconds.size(); ++op)
                    if (milliseconds[op] > 0)
                        std::clog << "[Vulkan profile] first_batch op=" << op
                                  << " gpu_ms=" << milliseconds[op] << '\n';
                p.profiled = true;
            }
            ++batches;
        }
        result.iterations += static_cast<int>(z[6]);
        for (uint32_t i = 0; i < p.n; ++i)
            x[i] += w[i];
        double residual2 = 0;
        for (uint32_t i = 0; i < p.n; ++i) {
            double r = solveRhs[i];
            for (uint32_t j = p.rows[i]; j < p.rows[i + 1]; ++j)
                r -= p.originalValues[j] * x[p.columns[j]];
            residual[i] = r;
            residual2 += r * r;
        }
        result.reduction = norm2 == 0 ? 0 : std::sqrt(residual2 / norm2);
        result.converged = std::isfinite(result.reduction) && result.reduction <= tolerance
            && std::all_of(x.begin(), x.end(), [](double v) { return std::isfinite(v); });
        if (result.converged || !std::isfinite(result.reduction) || z[6] == 0)
            break;
    }
    if (p.dilu) {
        std::vector<double> original(p.n);
        for (uint32_t i = 0; i < p.n / 3; ++i)
            for (uint32_t c = 0; c < 3; ++c)
                original[3 * p.permutation[i] + c] = x[3 * i + c];
        x.swap(original);
    }
    result.seconds = std::chrono::duration<double>(Clock::now() - start).count();
    return result;
}
} // namespace Reslab
