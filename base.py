import copy
class TreeNode:
    def __init__(self, name, children=None):
        self.name = name
        self.children = children or []
    def fx(self, fxname):
        return TreeNode("f_"+fxname, [self])
    def __repr__(self):
        return string_equation(str_form(self))
    def __eq__(self, other):
        return str_form(self)==str_form(other)
    def __add__(self, other):
        return TreeNode("f_add", [self,other])
    def __mul__(self, other):
        return TreeNode("f_mul", [self,other])
    def __sub__(self, other):
        return self + other * tree_form("d_-1")
    def __pow__(self, other):
        return TreeNode("f_pow", [self,other])
    def __truediv__(self, other):
        return self * other ** tree_form("d_-1")
    def __and__(self, other):
        return TreeNode("f_and", [self,other])
    def __or__(self, other):
        return TreeNode("f_or", [self,other])
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
    if equation_tree.name == "f_list":
        return "["+",".join([string_equation_helper(child) for child in equation_tree.children])+"]"
    s = "(" 
    if len(equation_tree.children) == 1 or equation_tree.name in ["f_integrate"]:
        s = equation_tree.name[2:] + s
    sign = {"f_cosec":"?" , "f_sec":"?", "f_cot": "?", "f_circumcenter":"?", "f_transpose":"?", "f_exp":"?", "f_abs":"?", "f_log":"?", "f_and":"&", "f_or":"|", "f_sub":"-", "f_neg":"?", "f_inv":"?", "f_add": "+", "f_mul": "*", "f_pow": "^", "f_poly": ",", "f_div": "/", "f_int": ",", "f_sub": "-", "f_dif": "?", "f_sin": "?", "f_cos": "?", "f_tan": "?", "f_eq": "=", "f_sqt": "?"}
    for child in equation_tree.children:
        k = None
        if equation_tree.name not in sign.keys():
            k = ","
        else:
            k = sign[equation_tree.name]
        s+= string_equation_helper(copy.deepcopy(child)) + k
    s = s[:-1] + ")"
    return s
def string_equation(eq):
    
    alpha = ["x", "y", "z"]+[chr(x+ord("a")) for x in range(0,23)]
    eq = tree_form(eq)
    
    for i, letter in enumerate(alpha):
        eq = replace(eq, tree_form("v_"+str(i)), tree_form(letter))
    eq = str_form(eq)
    eq = eq.replace("d_", "")
    eq = eq.replace("s_", "")
    eq = eq.replace("'", "")
    
    return string_equation_helper(tree_form(eq))
