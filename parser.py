import copy
from lark import Lark, Tree
from base import *
import re
def extract_nested_lists(s, n=0):
    stack = []
    results = []
    current = ''
    indices = []
    start = None
    for i, char in enumerate(s):
        if char == '[':
            if not stack:
                start = i
            else:
                current += char
            stack.append('[')
        elif char == ']':
            stack.pop()
            if stack:
                current += char
            else:
                results.append(current)
                indices.append((start, i + 1))
                current = ''
        elif stack:
            current += char
    formatted_results = [re.sub(r'([^\[\],\s]+)', r'"\1"', lst) for lst in results]
    nested_lists = [eval(f'[{lst}]', {'__builtins__': {}}, {}) for lst in formatted_results]
    for idx, (start, end) in enumerate(reversed(indices)):
        s = s[:start] + f'{chr(len(indices)-idx-1+ord("A")+n)}' + s[end:]
    return nested_lists, s
grammar = """
?start: expr
?expr: logic_or

?logic_or: logic_and
         | logic_or "|" logic_and  -> or

?logic_and: equality
          | logic_and "&" equality  -> and

?equality: arithmetic
         | equality "=" arithmetic  -> eq

?arithmetic: term
           | arithmetic "+" term   -> add
           | arithmetic "-" term   -> sub

?term: factor
     | term "*" factor  -> mul
     | term "/" factor  -> div

?factor: not_expr
       | factor "^" not_expr  -> pow

?not_expr: "!" not_expr  -> not
         | "-" not_expr  -> neg   // Handles negative numbers and expressions
         | base

?base: NUMBER            -> number
     | FUNC_NAME "(" expr ")" -> func
     | VARIABLE          -> variable
     | "(" expr ")"      -> paren

FUNC_NAME: "sin" | "circumcenter" | "cos" | "tan" | "log" | "sqrt" | "int" | "dif" | "abs" | "transpose" | "exp" | "cosec" | "sec" | "cot"

VARIABLE: "x" | "y" | "z" | "A" | "B" | "C" | "D" | "a" | "b" | "c" | "d" | "e" | "f"

%import common.NUMBER
%import common.WS_INLINE
%ignore WS_INLINE
"""


parser = Lark(grammar, start='start', parser='lalr')
def take_input(equation):
  equation = equation.replace(" ", "")
  mdata, equation = extract_nested_lists(equation)
  def conv_list(eq):
    if isinstance(eq, list):
        return TreeNode("f_list", [conv_list(child) for child in eq])
    else:
        return tree_form(take_input(eq))
  mdata = [conv_list(element) for element in mdata]
  parse_tree = parser.parse(equation)
  def convert_to_treenode(parse_tree):
      def tree_to_treenode(tree):
          if isinstance(tree, Tree):
              node = TreeNode(tree.data)
              node.children = [tree_to_treenode(child) for child in tree.children]
              return node
          else:
              return TreeNode(str(tree))
      return tree_to_treenode(parse_tree)
  def remove_past(equation):
      if equation.name in {"number", "paren", "func", "variable"}:
          if len(equation.children) == 1:
            for index, child in enumerate(equation.children):
              equation.children[index] = remove_past(child)
            return equation.children[0]
          else:
            for index, child in enumerate(equation.children):
              equation.children[index] = remove_past(child)
            return TreeNode(equation.children[0].name, equation.children[1:])
      coll = TreeNode(equation.name, [])
      for child in equation.children:
          coll.children.append(remove_past(child))
      return coll
  tree_node = convert_to_treenode(parse_tree)
  tree_node = remove_past(tree_node)
  def fxchange(tree_node):
    return TreeNode("f_"+tree_node.name if tree_node.name in ["cosec", "sec", "cot", "or", "not", "and", "exp", "circumcenter", "transpose", "eq", "sub", "neg", "inv", "add", "sin", "cos", "tan", "mul", "int", "dif", "pow", "div", "log", "abs"]\
                    else "d_"+tree_node.name, [fxchange(child) for child in tree_node.children])
  tree_node = fxchange(tree_node)
  for i in range(26):
    alpha = ["x", "y", "z"]+[chr(x+ord("a")) for x in range(0,23)]
    tree_node = replace(tree_node, tree_form("d_"+alpha[i]), tree_form("v_"+str(i)))
  for i in range(len(mdata)):
    tree_node = replace(tree_node, tree_form("d_"+chr(i+ord("A"))), mdata[i])
  tree_node = str_form(tree_node)
  
  return tree_node
