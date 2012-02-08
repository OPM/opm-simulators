import argparse
from xml.dom import minidom

# fuzzy compare XML tree from XML strings
def isFuzzyEqualXml(xml1, xml2, absolute, relative):
    dom1 = minidom.parseString(xml1)
    dom2 = minidom.parseString(xml2)
    return isFuzzyEqualNode(dom1.documentElement, dom2.documentElement, absolute, relative)

# fuzzy compare of XML nodes
def isFuzzyEqualNode(node1, node2, absolute, relative):
    if node1.tagName != node2.tagName:
        print 'The name of the node differs in ', node1.tagName, ' and ', node2.tagName
        return False
    if sorted(node1.attributes.items()) != sorted(node2.attributes.items()):
        print 'Attributs differ in node ', node1.tagName
        return False
    if len(node1.childNodes) != len(node2.childNodes):
        print 'Number of children differs in node ', node1.tagName
        return False
    for node1child, node2child in zip(node1.childNodes, node2.childNodes):
        if node1child.nodeType != node2child.nodeType:
            print 'Node type differs in ', note1.tagName
            return False
        if node1child.nodeType == node1child.TEXT_NODE and not isFuzzyEqualText(node1child.data, node2child.data, absolute, relative):
            print 'Data differs in node ', node1.tagName
            return False
        if node1child.nodeType == node1child.ELEMENT_NODE and not isFuzzyEqualNode(node1child, node2child, absolute, relative):
            return False
    return True

# fuzzy compare of text consisting of whitespace separated numbers
def isFuzzyEqualText(text1, text2, absolute, relative):
    list1 = text1.split()
    list2 = text2.split()
    # difference only in whitespace?
    if (list1 == list2):
        return True
    # compare number by number
    for number1, number2 in zip(list1, list2):
        number1 = float(number1)
        number2 = float(number2)
        if (abs(number1 - number2) > absolute 
            and (number2 == 0.0 or abs(abs(number1 / number2) - 1.0) > relative)):
            print 'Difference to large between', number1, ' and ', number2
            return False
    return True

# main programm
# handle arguments and print help message
parser = argparse.ArgumentParser(description='Fuzzy compare two VTK\
    (Visualization Toolkit) files. The files are accepted if for every\
    value the difference is below the absolute error or below the\
    relative error or below both.')
parser.add_argument('vtu_file_1', type=open,
    help='first file to compare')
parser.add_argument('vtu_file_2', type=open,
    help='second file to compare')
parser.add_argument('-r', '--relative', type=float, default=1e-2,
    help='maximum relative error (default=1e-2)')
parser.add_argument('-a', '--absolute', type=float, default=1e-9,
    help='maximum relative error (default=1e-9)')
args = parser.parse_args()

# fuzzy compare
if (isFuzzyEqualXml(args.vtu_file_1.read(), args.vtu_file_2.read(), args.absolute, args.relative)):
    exit
else:
    exit(1)


