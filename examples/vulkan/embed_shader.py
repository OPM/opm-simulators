#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
import pathlib
import struct
import sys

source, target = map(pathlib.Path, sys.argv[1:])
data = source.read_bytes()
words = struct.unpack("<" + "I" * (len(data) // 4), data)
target.write_text(
    f"#pragma once\n#include <cstdint>\ninline constexpr std::uint32_t {target.stem}[] = {{\n"
    + ",\n".join(",".join(hex(w) for w in words[i : i + 8]) for i in range(0, len(words), 8))
    + "\n};\n"
)
