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
    
    # Ensure all elements inside brackets are treated as strings
    formatted_results = [re.sub(r'([^\[\],\s]+)', r'"\1"', lst) for lst in results]
    
    # Convert extracted strings to Python lists safely
    nested_lists = [eval(f'[{lst}]', {'__builtins__': {}}, {}) for lst in formatted_results]
    
    # Replace original substrings with "v_n"
    for idx, (start, end) in enumerate(reversed(indices)):
        s = s[:start] + f'{chr(len(indices)-idx-1+ord("A")+n)}' + s[end:]
    
    return nested_lists, s

grammar = """
?start: expr

?expr: term
     | expr "+" term   -> add
     | expr "-" term   -> sub

?term: factor
     | term "*" factor  -> mul
     | term "/" factor  -> div

?factor: base
       | factor "^" base  -> pow

?base: NUMBER            -> number
     | FUNC_NAME "(" expr ")" -> func
     | VARIABLE          -> variable
     | "(" expr ")"      -> paren

FUNC_NAME: "sin" | "cos" | "tan" | "log" | "sqrt" | "int" | "dif"
VARIABLE: "x" | "y" | "z" | "A" | "B" | "C"

%import common.NUMBER
%import common.WS_INLINE
%ignore WS_INLINE
"""
parser = Lark(grammar, start='start', parser='lalr')
def take_input(equation):
  equation = equation.replace(" ", "")
  mdata, equation = extract_nested_lists(equation)
  
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
    return TreeNode("f_"+tree_node.name if tree_node.name in ["sub", "neg", "inv", "add", "sin", "cos", "tan", "mul", "int", "dif", "pow", "div"]\
                    else "d_"+tree_node.name, [fxchange(child) for child in tree_node.children])
  tree_node = fxchange(tree_node)
  
  for i in range(3):
    tree_node = replace(tree_node, tree_form("d_"+["x", "y", "z"][i]), tree_form("v_"+str(i)))
  for i in range(len(mdata)):
    tree_node = replace(tree_node, tree_form("d_"+chr(i+ord("A"))), tree_form("v_"+str(mdata[i]).replace(" ", "")))
  
  tree_node = str_form(tree_node)
  
  return tree_node
