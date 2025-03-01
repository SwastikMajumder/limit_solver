import copy
import itertools
from base import *

def flatten_tree(node):
    if not node.children:
        return node
    if node.name in ("f_add", "f_mul"):
        merged_children = []
        for child in node.children:
            flattened_child = flatten_tree(child)
            if flattened_child.name == node.name:
                merged_children.extend(flattened_child.children)
            else:
                merged_children.append(flattened_child)
        return TreeNode(node.name, merged_children)
    else:
        node.children = [flatten_tree(child) for child in node.children]
        return node
    
import parser
def convert_sub2neg(eq):
    def a1(eq):
        if eq.name == "f_sub":
            return TreeNode("f_add", [eq.children[0], TreeNode("f_neg", [eq.children[1]])])
        elif eq.name == "f_div":
            return TreeNode("f_mul", [eq.children[0], TreeNode("f_inv", [eq.children[1]])])
        term = TreeNode(eq.name, [])
        for child in eq.children:
            term.children.append(a1(child))
        return term

    def a2(eq):
        if eq.name == "f_neg":
            return TreeNode("f_mul", [tree_form("d_-1"), eq.children[0]])
        elif eq.name == "f_inv":
            return TreeNode("f_pow", [eq.children[0], tree_form("d_-1")])
        term = TreeNode(eq.name, [])
        for child in eq.children:
            term.children.append(a2(child))
        return term
    eq = tree_form(eq)
    while "f_sub" in str_form(eq) or "f_div" in str_form(eq):
        eq = a1(eq)
    eq = flatten_tree(eq)
    while "f_neg" in str_form(eq) or "f_inv" in str_form(eq):
        eq = a2(eq)
        
    return str_form(eq)

def simple(eq):
    if (eq.name[:2] == "f_" and eq.name != "f_add") or eq.name[:2] in ["d_","v_"]:
        eq = TreeNode("f_add", [eq])
    if eq.name == "f_add":
        dic = {}
        num3 = 0
        for i in range(len(eq.children)-1,-1,-1):
            if eq.children[i].name[:2] == "d_":
                num3 += int(eq.children[i].name[2:])
                eq.children.pop(i)
        for child in eq.children:
            num = 1
            if child.name == "f_mul":
                for i in range(len(child.children)-1,-1,-1):
                    if child.children[i].name[:2] == "d_":
                        num *= int(child.children[i].name[2:])
                        child.children.pop(i)
                dic2 = {}
                for i in range(len(child.children)-1,-1,-1):
                    num2 = 1
                    child2 = str_form(child.children[i])
                    if child.children[i].name == "f_pow" and child.children[i].children[1].name[:2]== "d_":
                        num2 = int(child.children[i].children[1].name[2:])
                        child2 = str_form(child.children[i].children[0])
                    if child2 not in dic2.keys():
                        dic2[child2] = num2
                    else:
                        dic2[child2] += num2
                newchild = TreeNode("f_mul", [])
                for key in sorted(dic2.keys()):
                    if dic2[key] == 0:
                        continue
                    if dic2[key] == 1:
                        newchild.children.append(tree_form(key))
                    else:
                        newchild.children.append(TreeNode("f_pow", [tree_form(key), tree_form("d_"+str(dic2[key]))]))
                if len(newchild.children)==1:
                    newchild = newchild.children[0]
                child = newchild
            if child.name in {"f_add", "f_mul"}:
                child = TreeNode(child.name, [tree_form(x) for x in sorted([str_form(x) for x in child.children])])
            if child.name == "f_mul" and child.children == []:
                num3 += num
            else:
                child = str_form(child)
                if child not in dic.keys():
                    dic[child] = num
                else:
                    dic[child] += num
        newchild = TreeNode("f_add", [])
        for key in sorted(dic.keys()):
            if dic[key] == 0:
                continue
            if dic[key] == 1:
                newchild.children.append(tree_form(key))
            else:
                newchild.children.append(TreeNode("f_mul", [tree_form("d_"+str(dic[key])), tree_form(key)]))
        for i in range(len(newchild.children)-1,-1,-1):
            if newchild.children[i].name[:2] == "d_":
                num3 += int(newchild.children[i].name[2:])
                newchild.children.pop(i)
        if num3 != 0:
            newchild.children.append(tree_form("d_"+str(num3)))
        if len(newchild.children)==1:
            newchild = newchild.children[0]
        elif len(newchild.children)==0:
            return tree_form("d_0")
        new = newchild
        if new.name in {"f_add", "f_mul"}:
            new = TreeNode(new.name, [tree_form(x) for x in sorted([str_form(x) for x in new.children])])
        new = flatten_tree(new)
        for i in range(len(new.children)):
            new.children[i] = simple(copy.deepcopy(new.children[i]))
        return new
def expand_eq(eq):
    eq = tree_form(eq)
    eq = simple(eq)
    return eq


def solve(eq):
    eq = convert_sub2neg(eq)
    eq= expand_eq(eq)
    eq = inversehandle(eq)
    return str_form(eq)
class Eq:
    def __init__(self, string, format_type="stringequation"):
        if format_type=="stringequation":
            self.equation= solve(parser.take_input(string))
        elif format_type=="treeform":
            self.equation= solve(str_form(string))
        elif format_type == "strform":
            self.equation= solve(string)
    def __eq__(self, other):
        x = (self-other).equation
        
        return x=="d_0"
    def __sub__(self, other):
        return Eq(solve(str_form(TreeNode("f_sub", [tree_form(self.equation),tree_form(other.equation)]))), "strform")
    def __add__(self, other):
        return Eq(solve(str_form(TreeNode("f_add", [tree_form(self.equation),tree_form(other.equation)]))), "strform")
    def __truediv__(self, other):
        return Eq(solve(str_form(TreeNode("f_div", [tree_form(self.equation),tree_form(other.equation)]))), "strform")
    def __mul__(self, other):
        return Eq(solve(str_form(TreeNode("f_mul", [tree_form(self.equation),tree_form(other.equation)]))), "strform")
    def factor(self):
        output = []
        eq = tree_form(self.equation)
        if eq.name != "f_mul":
            eq = TreeNode("f_mul", [eq])
        if eq.name == "f_mul":
            for child in eq.children:
                if child.name == "f_pow":
                    if int(child.children[1].name[2:]) == -1:
                        output.append(Eq("1")/Eq(child.children[0], "treeform"))
                    else:
                        for i in range(int(child.children[1].name[2:])):
                            output.append(Eq(child.children[0], "treeform"))
                else:
                    output.append(Eq(child, "treeform"))
        return output
    def __repr__(self):
        return string_equation(self.equation)

def substitute_val(eq, val):
    return Eq(replace(tree_form(eq.equation), tree_form("v_0"), tree_form("d_"+str(val))), "treeform")

def inversehandle(eq):
    if eq.name == "f_pow" and eq.children[0].name == "f_mul":
        arr = TreeNode("f_mul", [])
        for child in eq.children[0].children:
            arr.children.append(TreeNode("f_pow", [child, eq.children[1]]))
        return arr
    arr = TreeNode(eq.name, [])
    for child in eq.children:
        arr.children.append(inversehandle(child))
    return arr

def simplify(eq):
    if all(x.name[:2] == "d_" for x in eq.children):
        if eq.name == "f_pow":
            if int(eq.children[0].name[2:]) == 0 and int(eq.children[1].name[2:]) == 0:
                    return None
            if int(eq.children[1].name[2:]) > 0:
                return int(eq.children[0].name[2:]) ** int(eq.children[1].name[2:])
            if int(eq.children[1].name[2:]) == 0:
                return str(int(eq.children[0].name[2:])/abs(int(eq.children[0].name[2:])))
        elif eq.name == "f_add":
            arr = [int(x.name[2:]) for x in eq.children]
            return sum(arr)
        elif eq.name == "f_mul":
            arr = [int(x.name[2:]) for x in eq.children]
            p = 1
            for item in arr:
                p = p * item
            return p
    arr = TreeNode(eq.name, [])
    for child in eq.children:
        child = simplify(child)
        if child is None:
            return None
        if isinstance(child, int):
            child = tree_form("d_"+str(child))
        arr.children.append(child)
    return arr
    
                    
equation = None

while True:
    
    tmp = input(">>> ")
    try:
        orig = equation
        if tmp == "factor":
            for item in equation.factor():
                print(item)
        elif tmp == "show":
            print(equation)
        elif tmp == "limit 0":
            orig = equation
            con = True
            while con:
                con = False
                for item in equation.factor():
                    if tree_form(item.equation).name == "f_sin":
                        equation = equation/item
                        equation = equation*Eq(tree_form(item.equation).children[0], "treeform")
                        con = True
                        break
                    elif tree_form(item.equation).name == "f_pow" and tree_form(item.equation).children[1].name == "d_-1" and tree_form(item.equation).children[0].name == "f_sin":
                        equation = equation*Eq(tree_form(item.equation).children[0], "treeform")
                        equation = equation/Eq(tree_form(item.equation).children[0].children[0], "treeform")
                        con = True
                        break
                    elif tree_form(item.equation).name == "f_cos":
                        equation = equation/item
                        x = Eq(tree_form(item.equation).children[0], "treeform")
                        equation = equation*(Eq("1")-x*x/Eq("2"))
                        con = True
                        break
                    elif tree_form(item.equation).name == "f_pow" and tree_form(item.equation).children[1].name == "d_-1" and tree_form(item.equation).children[0].name == "f_cos":
                        equation = equation*Eq(tree_form(item.equation).children[0], "treeform")
                        x = Eq(tree_form(item.equation).children[0].children[0], "treeform")
                        equation = equation/(Eq("1")-x*x/Eq("2"))
                        con = True
                        break
            e = simplify(tree_form(substitute_val(equation, 0).equation))
            equation = Eq(e, "treeform")
            print("limit x->0 " + str(orig) + " = " + str(equation))
            equation = orig
        else:
            equation = Eq(tmp)
            print(equation)
    except Exception as error:
        equation = orig
        print(error)
