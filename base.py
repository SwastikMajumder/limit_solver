import copy

class TreeNode:
    def __init__(self, name, children=None):
        self.name = name
        self.children = children or []
    def __repr__(self):
        return string_equation(str_form(self))
def replace(equation, find, r):
  if str_form(equation) == str_form(find):
    return r
  col = TreeNode(equation.name, [])
  for child in equation.children:
    col.children.append(replace(child, find, r))
  return col
def tree_form(tabbed_strings):
    lines = tabbed_strings.split("\n")
    root = TreeNode("Root")
    current_level_nodes = {0: root}
    stack = [root]
    for line in lines:
        level = line.count(' ')
        node_name = line.strip()
        node = TreeNode(node_name)
        while len(stack) > level + 1:
            stack.pop()
        parent_node = stack[-1]
        parent_node.children.append(node)
        current_level_nodes[level] = node
        stack.append(node)
    return root.children[0]
def str_form(node):
    def recursive_str(node, depth=0):
        result = "{}{}".format(' ' * depth, node.name)
        for child in node.children:
            result += "\n" + recursive_str(child, depth + 1)
        return result
    return recursive_str(node)
def string_equation_helper(equation_tree):
    if equation_tree.children == []:
        return equation_tree.name
    s = "(" # bracket
    if len(equation_tree.children) == 1:
        s = equation_tree.name[2:] + s
    sign = {"f_sub":"-", "f_neg":"?", "f_inv":"?", "f_add": "+", "f_mul": "*", "f_pow": "^", "f_poly": ",", "f_div": "/", "f_int": ",", "f_sub": "-", "f_dif": "?", "f_sin": "?", "f_cos": "?", "f_tan": "?", "f_eq": "=", "f_sqt": "?"} # operation symbols
    for child in equation_tree.children:
        s+= string_equation_helper(copy.deepcopy(child)) + sign[equation_tree.name]
    s = s[:-1] + ")"
    return s
def string_equation(eq): 
    eq = eq.replace("v_0", "x")
    eq = eq.replace("v_1", "y")
    eq = eq.replace("v_2", "z")
    eq = eq.replace("d_", "")
    
    return string_equation_helper(tree_form(eq))
