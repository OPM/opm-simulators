import json
from sphinx.util.nodes import nested_parse_with_titles
from docutils.statemachine import ViewList
from sphinx.util.docutils import SphinxDirective
from docutils import nodes

def readDocStrings(directive, docstrings_path):
    print(docstrings_path)
    with open(docstrings_path, 'r') as file:
        docstrings = json.load(file)
    result = []
    for func_name, doc_info in docstrings.items():
        signature = doc_info.get('signature', '')
        docstring = doc_info.get('doc', '')

        # Create a ViewList instance for the function signature and docstring
        rst = ViewList()

        # Check if signature exists and prepend it to the docstring
        signature_line = f".. py:function:: {signature}" if signature else f".. py:function:: {func_name}()"
        rst.append(signature_line, source="")
        rst.append("", source="")
        # Add the docstring text if it exists
        if docstring:
            for line in docstring.split('\n'):
                rst.append(f"   {line}", source="")

        # Create a node that will be populated by nested_parse_with_titles
        node = nodes.section()
        node.document = directive.state.document
        # Parse the rst content
        nested_parse_with_titles(directive.state, rst, node)

        result.extend(node.children)
    return result

class SimulatorsDirective(SphinxDirective):
    required_arguments = 0
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec = {}

    def run(self):
        return readDocStrings(self, self.state.document.settings.env.app.config.opm_simulators_docstrings_path)

class CommonDirective(SphinxDirective):
    required_arguments = 0
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec = {}

    def run(self):
        return readDocStrings(self, self.state.document.settings.env.app.config.opm_common_docstrings_path)

def setup(app):
    app.add_config_value('opm_simulators_docstrings_path', None, 'env')
    app.add_config_value('opm_common_docstrings_path', None, 'env')
    app.add_directive("opm_simulators_docstrings", SimulatorsDirective)
    app.add_directive("opm_common_docstrings", CommonDirective)
