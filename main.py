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
    #eq = tree_form(eq)
    eq = simple(eq)
    return eq

def common(eq):
    
    if eq.name == "f_add":
        con = []
        for child in eq.children:
            if child.name == "f_mul":
                num = []
                den = []
                for child2 in child.children:
                    if child2.name == "f_pow" and child2.children[1].name[:2] == "d_" and int(child2.children[1].name[2:])<0:
                        n = int(child2.children[1].name[2:])
                        if n == -1:
                            den.append(child2.children[0])
                        else:
                            den.append(TreeNode("f_pow", [child2.children[0], tree_form("d_"+str(-n))]))
                    else:
                        num.append(child2)
                con.append([num, den])
            else:
                con.append([[child], []])
        if len(con)>1 and any(x[1] != [] for x in con):
            a = TreeNode("f_add", [])
            for i in range(len(con)):
                b = TreeNode("f_mul", [])
                if con[i][0] != []:
                    b.children += con[i][0]
                for j in range(len(con)):
                    if i ==j:
                        continue
                    b.children +=  con[j][1]
                
                if len(b.children) == 1:
                    a.children.append(b.children[0])
                elif len(b.children) > 1:
                    a.children.append(b)
            c = TreeNode("f_mul", [])
            for i in range(len(con)):
                c.children += con[i][1]
            if len(c.children)==1:
                c = c.children[0]
            c = TreeNode("f_pow", [c, tree_form("d_-1")])
            return TreeNode("f_mul", [a,c])
    arr = TreeNode(eq.name, [])
    for child in eq.children:
        arr.children.append(common(child))
    return arr
                        
def solve(eq):
    eq = convert_sub2neg(eq)
    
    eq = common(tree_form(eq))
    
    eq= expand_eq(eq)
    eq = inversehandle(eq)
    eq = inversehandle(eq)
    eq = flatten_tree(eq)
    eq= expand_eq(eq)
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
    if eq.name == "f_pow":
        base, exponent = eq.children
        if base.name == "f_pow":
            n = int(exponent.name[2:]) * int(base.children[1].name[2:])
            return base.children[0] if n == 1 else TreeNode("f_pow", [base.children[0], tree_form(f"d_{n}")])
        elif base.name == "f_mul":
            return TreeNode("f_mul", [TreeNode("f_pow", [child, exponent]) for child in base.children])
    arr = TreeNode(eq.name, [inversehandle(child) for child in eq.children])
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
            if abs(int(eq.children[0].name[2:])) == 1 and int(eq.children[1].name[2:]) == -1:
                return eq.children[0]
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

def replace_eq(eq):
    if eq.name=="f_tan":
        return TreeNode("f_div", [TreeNode("f_sin", [copy.deepcopy(eq.children[0])]), TreeNode("f_cos", [copy.deepcopy(eq.children[0])])])
    arr = TreeNode(eq.name, [])
    for child in eq.children:
        arr.children.append(replace_eq(child))
    return arr

def take_common(eq):
    lst = None
    orig = eq
    eq = tree_form(eq.equation)
    if eq.name == "f_mul":
        eq = [x for x in eq.children if x.name == "f_add"]
    elif eq.name == "f_add":
        eq = [eq]
    else:
        return orig
    eqlst = copy.deepcopy(eq)
    output = Eq("1")
    final = Eq("1")
    if tree_form(orig.equation).name == "f_mul":
        for item in [x for x in tree_form(orig.equation).children if x.name != "f_add"]:
            final = final * Eq(item, "treeform")
    for eq in eqlst:
        for child in eq.children:
            child = Eq(child, "treeform")
            if lst is None:
                lst = child.factor()
            else:
                lst2 = child.factor()
                for i in range(len(lst)-1,-1,-1):
                    if lst[i] not in lst2:
                        lst.pop(i)
                    else:
                        for j in range(len(lst2)):
                            if lst[i] == lst2[j]:
                                lst2.pop(j)
                                break
        ans = Eq("1")
        for item in lst:
            ans = ans * item
        final2 = Eq("0")
        for child in eq.children:
            final2 += Eq(child, "treeform")/ans
        final = final2 * final
        output = output * ans
    output = output*final
    return Eq(simplify(tree_form(output.equation)), "treeform")    
        
equation = None

def expand(eq):
    if eq.name == "f_mul":
        addchild = [[Eq(child2, "treeform") for child2 in child.children] for child in eq.children if child.name == "f_add"]
        otherchild = [child for child in eq.children if child.name != "f_add"]
        if len(otherchild) == 1 and otherchild[0].name[:2] == "d_" and len(addchild) == 1 and isinstance(addchild, list):
            add= Eq("0")
            for item in addchild[0]:
                add += Eq(otherchild[0], "treeform")*item
            return tree_form(add.equation)

    return TreeNode(eq.name, [expand(child) for child in eq.children])

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

            def approx(eq):
                if eq.name == "f_sin":
                    return eq.children[0]
                elif eq.name == "f_cos":
                    x = Eq(eq.children[0], "treeform")
                    return tree_form((Eq("1") - (x * x / Eq("2"))).equation)
                arr = TreeNode(eq.name, [])
                for child in eq.children:
                    arr.children.append(approx(child))
                return arr
            equation = Eq(approx(tree_form(equation.equation)), "treeform")
            equation = Eq(approx(tree_form(equation.equation)), "treeform")
            equation = Eq(expand(tree_form(equation.equation)), "treeform")
            equation = Eq(expand(tree_form(equation.equation)), "treeform")
            equation2 = equation
            e = simplify(tree_form(substitute_val(equation, 0).equation))
            if isinstance(e, int):
                e = tree_form("d_"+str(e))
            equation = Eq(e, "treeform")
            print("limit x->0 " + str(orig) + " = " + str(equation2) + " = " + str(equation))
            equation = orig
        else:
            equation = Eq(tmp)
            equation = Eq(replace_eq(tree_form(equation.equation)), "treeform")
            
            equation = take_common(equation)
            print(equation)
    except Exception as error:
        equation = orig
        print(error)
