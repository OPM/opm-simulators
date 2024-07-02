import json
from sphinx.util.nodes import nested_parse_with_titles
from docutils.statemachine import ViewList
from sphinx.util.docutils import SphinxDirective
from docutils import nodes

def read_doc_strings(directive, docstrings_path):
    print(docstrings_path)
    with open(docstrings_path, 'r') as file:
        docstrings = json.load(file)
    result = []
    for name, item in docstrings.items():
        # Create a ViewList instance for the function signature and docstring
        rst = ViewList()

        # Check if signature exists and prepend it to the docstring
        signature = item.get('signature', '')
        item_type = item.get('type', 'method')
        signature_line = f".. py:{item_type}:: {signature}" if signature else f".. py:{item_type}:: {name}()"
        rst.append(signature_line, source="")
        rst.append("", source="")

        # Add the docstring text if it exists
        docstring = item.get('doc', '')
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
        return read_doc_strings(self, self.state.document.settings.env.app.config.opm_simulators_docstrings_path)

class CommonDirective(SphinxDirective):
    required_arguments = 0
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec = {}

    def run(self):
        return read_doc_strings(self, self.state.document.settings.env.app.config.opm_common_docstrings_path)

def setup(app):
    app.add_config_value('opm_simulators_docstrings_path', None, 'env')
    app.add_config_value('opm_common_docstrings_path', None, 'env')
    app.add_directive("opm_simulators_docstrings", SimulatorsDirective)
    app.add_directive("opm_common_docstrings", CommonDirective)
