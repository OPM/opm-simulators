import json
import sys

def generate_hpp_from_json(json_path: str, output_hpp_path: str):
    with open(json_path, 'r', encoding='utf-8') as file:
        docstrings = json.load(file)

    hpp_content = """\
#ifndef PYBLACKOILSIMULATORDOC_HPP
#define PYBLACKOILSIMULATORDOC_HPP

// Generated docstrings for PyBlackOilSimulator
namespace Opm::Pybind::DocStrings {
"""

    for func_name, info in docstrings.items():
        signature = info['signature']
        doc = info['doc'].replace('\n', '\n    ')
        hpp_content += f"""
static constexpr char {func_name}_docstring[] = R\"doc(
{doc}
)doc\";\n"""

    hpp_content += """\
} // namespace Opm::Pybind::DocStrings

#endif // PYBLACKOILSIMULATORDOC_HPP
"""

    with open(output_hpp_path, 'w', encoding='utf-8') as file:
        file.write(hpp_content)

if __name__ == "__main__":
    # Check that exactly two command line arguments are provided
    if len(sys.argv) != 3:
        print("Usage: python generate_docstring_hpp.py <json_path> <output_hpp_path>")
        sys.exit(1)
    # Extract json_path and output_hpp_path from command line arguments
    json_path = sys.argv[1]
    output_hpp_path = sys.argv[2]
    generate_hpp_from_json(json_path, output_hpp_path)
