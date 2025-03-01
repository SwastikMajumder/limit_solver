import copy
from lark import Lark, Tree
from base import *

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
VARIABLE: "x" | "y" | "z"

%import common.NUMBER
%import common.WS_INLINE
%ignore WS_INLINE
"""
parser = Lark(grammar, start='start', parser='lalr')
def take_input(equation):
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
  tree_node = str_form(tree_node)
  
  for item in ["sub", "neg", "inv", "add", "sin", "cos", "tan", "mul", "int", "dif", "pow", "div"]:
    tree_node = tree_node.replace(str(item), "f_" + str(item))
  tree_node = tree_form(tree_node)
  for i in range(100,-1,-1):
    tree_node = replace(tree_node, tree_form(str(i)), tree_form("d_"+str(i)))
  for i in range(3):
    tree_node = replace(tree_node, tree_form(["x", "y", "z"][i]), tree_form("v_"+str(i)))
  tree_node = str_form(tree_node)
  return tree_node
