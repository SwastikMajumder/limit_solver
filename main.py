import copy
import itertools
from base import *
from fractions import Fraction

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
                num3 += Fraction(eq.children[i].name[2:])
                eq.children.pop(i)
            elif eq.children[i].name == "f_pow" and eq.children[i].children[0].name[:2]=="d_" and Fraction(eq.children[i].children[1].name[2:]) < 0:
                num3 += 1/(Fraction(eq.children[i].children[0].name[2:])**(-Fraction(eq.children[i].children[1].name[2:])))
                eq.children.pop(i)
        for child in eq.children:
            num = 1
            if child.name == "f_mul":
                for i in range(len(child.children)-1,-1,-1):
                    if child.children[i].name[:2] == "d_":
                        num *= Fraction(child.children[i].name[2:])
                        child.children.pop(i)
                    elif child.children[i].name == "f_pow" and child.children[i].children[0].name[:2]=="d_" and Fraction(child.children[i].children[1].name[2:]) < 0:
                        num *= 1/(Fraction(child.children[i].children[0].name[2:])**(-Fraction(child.children[i].children[1].name[2:])))
                        child.children.pop(i)
                dic2 = {}
                for i in range(len(child.children)-1,-1,-1):
                    num2 = 1
                    child2 = str_form(child.children[i])
                    if child.children[i].name == "f_pow" and child.children[i].children[1].name[:2]== "d_":
                        num2 = Fraction(child.children[i].children[1].name[2:])
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
                num3 += Fraction(newchild.children[i].name[2:])
                newchild.children.pop(i)
        if num3 != 0:
            x = tree_form("d_"+str(num3))
            newchild.children.append(x)
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
    def frac(eq):
        if eq.name[:2] == "d_":
            n = Fraction(eq.name[2:])
            if n.denominator == 1:
                return tree_form("d_"+str(n.numerator))
            elif n.numerator == 1:
                b = tree_form("d_"+str(n.denominator))
                b = TreeNode("f_pow", [b, tree_form("d_-1")])
                return b
            else:
                a = tree_form("d_"+str(n.numerator))
                b = tree_form("d_"+str(n.denominator))
                b = TreeNode("f_pow", [b, tree_form("d_-1")])
                return TreeNode("f_mul", [a,b])
        return TreeNode(eq.name, [frac(child) for child in eq.children])
    return frac(eq)

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
                else:
                    a.children.append(tree_form("d_1"))
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
def common2(eq):
    eq = common(eq)
    return TreeNode(eq.name, [common2(child) for child in eq.children])
def solve(eq):
    eq = convert_sub2neg(eq)
    eq = tree_form(eq)
    eq = expand_eq(eq)
    eq = inversehandle(eq)
    eq = inversehandle(eq)
    eq = flatten_tree(eq)
    eq = simplify(eq)
    eq = expand_eq(eq)
    eq = simp(eq)
    return str_form(eq)

#1/tan(x)-tan(x)-2/tan((2*x))
def simp(eq):
    if eq.name[2:] in "sin cos tan".split(" "):
        
        eq2 = str_form(eq.children[0])
        eq2 = convert_sub2neg(eq2)
        eq2 = tree_form(eq2)
        eq2 = expand_eq(eq2)
        eq2 = flatten_tree(eq2)
        eq2 = simplify(eq2)
        eq2 = expand_eq(eq2)
        
        eq.children[0] = eq2
        return TreeNode(eq.name, [eq.children[0]])
    elif eq.children != []:
        return TreeNode(eq.name, [simp(child) for child in eq.children])
    else:
        return eq
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
                    if int(child.children[1].name[2:]) < 0:
                        for i in range(-int(child.children[1].name[2:])):
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
    eq = replace(tree_form(eq.equation), tree_form("v_0"), tree_form("d_"+str(val)))
    return eq

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
    if eq is None:
        return None
    if all(x.name[:2] == "d_" for x in eq.children):
        if eq.name == "f_sin":
            if eq.children[0].name == "d_0":
                return tree_form("d_0")
        if eq.name == "f_cos":
            if eq.children[0].name == "d_0":
                return tree_form("d_1")
        if eq.name == "f_pow":
            if int(eq.children[0].name[2:]) == 0 and int(eq.children[1].name[2:]) <= 0:
                return None
            if int(eq.children[1].name[2:]) > 0:
                n = int(eq.children[0].name[2:]) ** int(eq.children[1].name[2:])
                return tree_form("d_"+str(n))
            if int(eq.children[1].name[2:]) == 0:
                return None
            if abs(int(eq.children[0].name[2:])) == 1 and int(eq.children[1].name[2:]) == -1:
                return eq.children[0]
        elif eq.name == "f_add":
            arr = [int(x.name[2:]) for x in eq.children]
            return tree_form("d_"+str(sum(arr)))
        elif eq.name == "f_mul":
            arr = [int(x.name[2:]) for x in eq.children]
            p = 1
            for item in arr:
                p = p * item
            return tree_form("d_"+str(p))
    arr = TreeNode(eq.name, [])
    for child in eq.children:
        child = simplify(child)
        if child is None:
            return None
        arr.children.append(child)
    return arr

def replace_eq(eq):
    if eq.name[2:] in ["tan", "sin", "cos"]:
        eq = TreeNode(eq.name, [tree_form(Eq(eq.children[0], "treeform").equation)])
    
    if eq.name=="f_tan":
        return replace_eq(TreeNode("f_div", [TreeNode("f_sin", [copy.deepcopy(eq.children[0])]), TreeNode("f_cos", [copy.deepcopy(eq.children[0])])]))
    if eq.name == "f_cos":
        child = Eq(eq.children[0], "treeform")
        child2 = child*Eq("d_-1", "strform")
        if child.equation.count("d_-") > child2.equation.count("d_-"):
            eq = TreeNode("f_cos", [tree_form(child2.equation)])
    elif eq.name == "f_sin":
        child = Eq(eq.children[0], "treeform")
        child2 = child*Eq("d_-1", "strform")
        if child.equation.count("d_-") > child2.equation.count("d_-"):
            eq = TreeNode("f_mul", [tree_form("d_-1"), TreeNode("f_sin", [tree_form(child2.equation)])])
    if eq.name in ["f_sin", "f_cos"] and eq.children[0].name == "f_add" and len(eq.children[0].children)==2:
        if eq.name == "f_sin":
            
            a = TreeNode("f_sin", [eq.children[0].children[0]])
            b = TreeNode("f_cos", [eq.children[0].children[1]])
            c = TreeNode("f_cos", [eq.children[0].children[0]])
            d = TreeNode("f_sin", [eq.children[0].children[1]])
            a, b, c, d = [Eq(x, "treeform") for x in [a, b, c, d]]
            
            return replace_eq(tree_form((a*b+c*d).equation))
        
        elif eq.name == "f_cos":
            a = TreeNode("f_cos", [eq.children[0].children[0]])
            b = TreeNode("f_cos", [eq.children[0].children[1]])
            c = TreeNode("f_sin", [eq.children[0].children[0]])
            d = TreeNode("f_sin", [eq.children[0].children[1]])
            a, b, c, d = [Eq(x, "treeform") for x in [a, b, c, d]]
            return replace_eq(tree_form((a*b-c*d).equation))
    if eq.name in ["f_sin", "f_cos"] and diffany2(Eq(eq.children[0], "treeform")).equation == "d_2":
        child = tree_form((Eq(eq.children[0], "treeform")/Eq("2")).equation)
        a = TreeNode("f_sin", [child])
        b = TreeNode("f_cos", [child])
        a, b = [Eq(x, "treeform") for x in [a, b]]
        if eq.name == "f_sin":
            return tree_form((Eq("2")*a*b).equation)
        else:
            return tree_form((b*b-a*a).equation)
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
    return output 
        


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

def expand2(eq):
    if eq.name == "f_mul":
        addchild = [[Eq(child2, "treeform") for child2 in child.children] for child in eq.children if child.name == "f_add"]
        otherchild = [Eq(child, "treeform") for child in eq.children if child.name != "f_add"]
        other = Eq("1")
        for item in otherchild:
            other= other * item
        if len(addchild) == 1 and isinstance(addchild, list):
            add= Eq("0")
            for item in addchild[0]:
                add = add + other*item
            return tree_form(add.equation)

    return TreeNode(eq.name, [expand2(child) for child in eq.children])

def numdem(equation):
    num = Eq("1")
    den = Eq("1")
    for item in equation.factor():
        t = tree_form(item.equation)
        if t.name == "f_pow" and int(t.children[1].name[2:]) < 0:
            den = den*item
        else:
            num = num*item
    return [num, Eq("1")/den]

def diff(equation):
    eq = tree_form(equation.equation)
    if "v_" not in str_form(eq):
        return Eq("0")
    if eq.name == "f_add":
        add = Eq("0")
        for child in eq.children:
            add += diff(Eq(child, "treeform"))
        return add
    elif eq.name == "f_mul":
        add = Eq("0")
        for i in range(len(eq.children)):
            new = copy.deepcopy(eq)
            new.children.pop(i)
            add += diff(Eq(eq.children[i], "treeform"))*Eq(new, "treeform")
        return add
    elif eq.name == "f_sin":
        eq.name = "f_cos"
        return diff(Eq(eq.children[0], "treeform"))*Eq(eq, "treeform")
    elif eq.name == "f_cos":
        eq.name = "f_sin"
        return Eq("d_-1", "strform")*diff(Eq(eq.children[0], "treeform"))*Eq(eq, "treeform")
    elif eq.name[:2] == "v_":
        return Eq(TreeNode("f_dif", [eq]), "treeform")
    elif eq.name == "f_pow" and eq.children[1].name[:2] == "d_":
        base, power = eq.children
        dbase = diff(Eq(base, "treeform"))
        b1 = Eq(power, "treeform") - Eq("1")
        bab1 = TreeNode("f_pow", [base, tree_form(b1.equation)])
        return Eq(power, "treeform")*Eq(bab1, "treeform")*dbase
def diffx2(equation):
    eq = tree_form(equation.equation)
    if eq.name == "f_dif":
        if eq.children[0].name == "v_0":
            return Eq("1")
        return Eq("0")
    return Eq(TreeNode(eq.name, [tree_form(diffx2(Eq(child, "treeform")).equation) for child in eq.children]), "treeform")
def diffany(equation):
    eq = tree_form(equation.equation)
    if eq.name == "f_dif":
        if eq.children[0].name[:2] == "v_":
            return Eq("1")
    return Eq(TreeNode(eq.name, [tree_form(diffany(Eq(child, "treeform")).equation) for child in eq.children]), "treeform")
def diffx(equation):
    equation = diff(equation)
    equation = diffx2(equation)
    return equation
def diffany2(equation):
    equation = diff(equation)
    equation = diffany(equation)
    return equation
def approx(equation):
    con = True
    while con:
        con = False
        for item in equation.factor():
            if tree_form(item.equation).name == "f_sin":
                if diffx(diffx(Eq(tree_form(item.equation).children[0], "treeform"))) != Eq("0"):
                    continue
                
                equation = equation/item
                equation = equation*Eq(tree_form(item.equation).children[0], "treeform")
                con = True
                break
            elif tree_form(item.equation).name == "f_pow" and tree_form(item.equation).children[1].name == "d_-1" and tree_form(item.equation).children[0].name == "f_sin":
                if diffx(diffx(Eq(tree_form(item.equation).children[0].children[0], "treeform"))) != Eq("0"):
                    continue
                equation = equation*Eq(tree_form(item.equation).children[0], "treeform")
                equation = equation/Eq(tree_form(item.equation).children[0].children[0], "treeform")
                con = True
                break
            elif tree_form(item.equation).name == "f_cos":
                if diffx(diffx(Eq(tree_form(item.equation).children[0], "treeform"))) != Eq("0"):
                    continue
                equation = equation/item
                x = Eq(tree_form(item.equation).children[0], "treeform")
                equation = equation*(Eq("1")-x*x/Eq("2"))
                con = True
                break
            elif tree_form(item.equation).name == "f_pow" and tree_form(item.equation).children[1].name == "d_-1" and tree_form(item.equation).children[0].name == "f_cos":
                if diffx(diffx(Eq(tree_form(item.equation).children[0].children[0], "treeform"))) != Eq("0"):
                    continue
                equation = equation*Eq(tree_form(item.equation).children[0], "treeform")
                x = Eq(tree_form(item.equation).children[0].children[0], "treeform")
                equation = equation/(Eq("1")-x*x/Eq("2"))
                con = True
                break
    return equation
def subslimit(equation):
    equation = Eq(expand2(tree_form(equation.equation)), "treeform")
    equation = simplify(tree_form(equation.equation))
    if equation is None:
        return None
    equation = Eq(equation, "treeform")
    equation = substitute_val(equation, 0)
    equation = simplify(equation)
    equation = simplify(equation)
    if equation is None:
        return None
    return Eq(equation, "treeform")
def lhostpital(equation):
    equation = simplify(tree_form(equation.equation))
    
    if equation is None:
        return None
    equation = Eq(equation, "treeform")
    e = substitute_val(equation, 0)
    e = simplify(e)
    e = simplify(e)
    if e is None:
        n, d = numdem(equation)
        ans1 = subslimit(n)
        ans2 = subslimit(d)
        if ans1 is not None and ans2 is not None and ans1 == Eq("0") and ans2 == Eq("0"):
            equation = diffx(n)/diffx(d)
            return equation
        return None
    return None
def additionlimit(equation):
    equation = Eq(expand2(tree_form(equation.equation)), "treeform")

    final = Eq("0")
    equation = tree_form(equation.equation)
    if equation.name == "f_add":
        for child in equation.children:
            tmp = find_limit(Eq(child, "treeform"))
            if tmp is not None:
                final += tmp
            else:
                final = None
                break
    else:
        final = None
    return final

def formula_1(eq):
    if eq.name == "f_pow" and abs(int(eq.children[1].name[2:])) == 2 and eq.children[0].name == "f_cos" and diffany2(Eq(eq.children[0].children[0], "treeform")).equation == "d_1":
        x = Eq(TreeNode("f_sin", [eq.children[0].children[0]]), "treeform")
        if int(eq.children[1].name[2:]) == 2:
            x = Eq("1")-x*x
        else:
            x = Eq("1")/(Eq("1")-x*x)
        return tree_form(x.equation)
    elif eq.name == "f_pow" and abs(int(eq.children[1].name[2:])) == 3 and eq.children[0].name == "f_cos" and diffany2(Eq(eq.children[0].children[0], "treeform")).equation == "d_1":
        x = Eq(TreeNode("f_sin", [eq.children[0].children[0]]), "treeform")
        y = Eq(TreeNode("f_cos", [eq.children[0].children[0]]), "treeform")
        if int(eq.children[1].name[2:]) == 3:
            x = y-x*x*y
        else:
            x = Eq("1")/(y-x*x*y)
        return tree_form(x.equation)
    return TreeNode(eq.name, [formula_1(child) for child in eq.children])
def approx_limit(equation):
    orig = equation
    equation = approx(equation)
    if orig == equation:
        return None
    return equation
def find_limit(equation,depth=3):
    ans = subslimit(equation)
    if ans is not None and "v_0" not in ans.equation:
        if "v_" in equation.equation.replace("v_0", ""):
            if "v_0" not in equation.equation:
                return equation
        else:
            return ans
    if depth == 0:
        return None
    
    for item in [lhostpital,additionlimit,approx_limit]:
        ans = item(equation)
        if ans is not None:
            out = find_limit(ans, depth-1)
            if out is not None:
                return out
    return None

equation = None
while True:
    tmp = input(">>> ")
    try:
        orig = equation
        if tmp == "factor":
            for item in equation.factor():
                print(item)
        elif tmp == "expand":
             #equation = Eq(flatten_tree(common(tree_form(equation.equation))), "treeform")
             equation = Eq(expand2(tree_form(equation.equation)), "treeform")
             print(equation)
        elif tmp == "show":
            print(equation)
        elif tmp == "fraction":
            equation = Eq(common2(tree_form(equation.equation)), "treeform")
            equation = Eq(expand(tree_form(equation.equation)), "treeform")
            equation = Eq(formula_1(tree_form(equation.equation)), "treeform")
            print(equation)
        elif tmp == "d/dx":
            print(diffx(equation))
        elif tmp == "limit 0":
            
            print("limit x->0 " + str(equation) + " = " + str(find_limit(equation)))
        else:
            equation = Eq(tmp,"stringequation")
            equation = Eq(replace_eq(tree_form(equation.equation)), "treeform")
            
            equation = take_common(equation)
            print(equation)
    except Exception as error:
        equation = orig
        print(error)
