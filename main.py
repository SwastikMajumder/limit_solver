import math
import copy
import itertools
from base import *
from fractions import Fraction
import parser
def perfect_square_root(n):
    if n < 0:
        root = -math.isqrt(-n)
    else:
        root = math.isqrt(n)
    return root if root * root == abs(n) else None
def pf(n):
    if not isinstance(n, Fraction):
        n = Fraction(n)
    a = n.numerator
    b = n.denominator
    if perfect_square_root(a) is not None and perfect_square_root(b) is not None:
        return Fraction(perfect_square_root(a),perfect_square_root(b))
    return None
def int2(string):
    tmp = Fraction(string)
    if tmp.denominator == 1:
        return tmp.numerator
    return tmp
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
def is_str_n(s, eq=None):
    s = tree_form(s)
    if s.name[:2] == "d_":
        if eq is None:
            return int2(s.name[2:])
        else:
            return int2(s.name[2:])==eq
    if eq is None:
        return None
    return False
def calc(eq):
    if eq.name[:2] == "d_":
        return Fraction(eq.name[2:])
    elif eq.name == "f_pow":
        a = calc(eq.children[0])
        b = calc(eq.children[1])
        if a is None or b is None:
            return None
        if b==Fraction(1,2):
            return pf(a)
        if b==Fraction(-1,1):
            return Fraction(1,1)/a
        return None
    elif eq.name == "f_mul":
        p = 1
        for child in eq.children:
            tmp = calc(child)
            if tmp is not None:
                p *= tmp
            else:
                return None
        return p
    elif eq.name == "f_add":
        p = 0
        for child in eq.children:
            
            tmp  = calc(child)
            
            if tmp is not None:
                
                p += tmp
            else:
                return None
        return p
    return None

def simple(eq):
    if (eq.name[:2] == "f_" and eq.name != "f_add") or eq.name[:2] in ["d_","v_"]:
        eq = copy.deepcopy(TreeNode("f_add", [eq]))
    
    if eq.name == "f_add":
        dic = {}
        num3 = Fraction(0)
        
        for i in range(len(eq.children)-1,-1,-1):
            tmp = calc(eq.children[i])
            if tmp is None:
                continue
            num3 += tmp
            eq.children.pop(i)
            
        for child in eq.children:
            num = 1
            if child.name == "f_mul":
                for i in range(len(child.children)-1,-1,-1):
                    tmp = calc(child.children[i])
                    if tmp is None:
                        continue
                    num *= tmp
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
    
def frac3(n):
    if n.denominator == 1:
        return tree_form("d_"+str(n.numerator))
    elif abs(n.numerator) == 1:
        b = tree_form("d_"+str(n.denominator*n.numerator))
        return TreeNode("f_pow", [b, tree_form("d_-1")])
    else:
        a = tree_form("d_"+str(n.numerator))
        b = tree_form("d_"+str(n.denominator))
        b = TreeNode("f_pow", [b, tree_form("d_-1")])
        return TreeNode("f_mul", [a,b])
def frac2(eq):
    if eq.name[:2] == "d_":
        return frac3(Fraction(eq.name[2:]))
    return TreeNode(eq.name, [frac2(child) for child in eq.children])
def expand_eq(eq):
    eq = simple(eq)
    eq = simple(eq)
    
    if "/" in str_form(eq):
        eq = frac2(eq)
        
        eq = simplify(eq)
        
        eq = simplify(eq)
    
    return eq

def common(eq):
    
    if eq.name == "f_add":
        con = []
        for child in eq.children:
            if child.name == "f_pow" and child.children[1].name[:2] == "d_" and int2(child.children[1].name[2:])<0:
                den = []
                n = int2(child.children[1].name[2:])
                if n == -1:
                    den.append(child.children[0])
                else:
                    den.append(TreeNode("f_pow", [child.children[0], tree_form("d_"+str(-n))]))
                con.append([[], den])
            elif child.name == "f_mul":
                num = []
                den = []
                for child2 in child.children:
                    if child2.name == "f_pow" and child2.children[1].name[:2] == "d_" and int2(child2.children[1].name[2:])<0:
                        n = int2(child2.children[1].name[2:])
                        
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
    #1/(x^8*(1+1/x^7)^(6/7))
    
    eq = convert_sub2neg(eq)
    eq = tree_form(eq)
    eq = expand_eq(eq)
    
    eq = flatten_tree(eq)
    
    eq = inversehandle(eq)
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
    @staticmethod
    def convert(term):
        return [Eq("d_" + str(item), "strform") if isinstance(item, int) else item for item in term]
    def __neg__(self):
        return Eq(solve(str_form(TreeNode("f_mul", [tree_form(self.equation),tree_form("d_-1")]))), "strform")
    def __sub__(self, other):
        
        return Eq(solve(str_form(TreeNode("f_sub", [tree_form(self.equation),tree_form(other.equation)]))), "strform")
    def __add__(self, other):
        return Eq(solve(str_form(TreeNode("f_add", [tree_form(self.equation),tree_form(other.equation)]))), "strform")
    def __truediv__(self, other):
        return Eq(solve(str_form(TreeNode("f_div", [tree_form(self.equation),tree_form(other.equation)]))), "strform")
    def __mul__(self, other):
        return Eq(solve(str_form(TreeNode("f_mul", [tree_form(self.equation),tree_form(other.equation)]))), "strform")
    def __pow__(self, other):
        self, other = Eq.convert([self, other])
        return Eq(solve(str_form(TreeNode("f_pow", [tree_form(self.equation),tree_form(other.equation)]))), "strform")
    def __rmul__(self, other):
        return self.__mul__(Eq("d_"+str(other), "strform"))
    def __radd__(self, other):
        return self.__add__(Eq("d_"+str(other), "strform"))
    def __rpow__(self, other):
        return self.__pow__(Eq("d_"+str(other), "strform"))
    def fx(self, fxname):
        return Eq(TreeNode("f_"+fxname, [tree_form(self.equation)]), "treeform")
    def factor(self):
        output = []
        eq = tree_form(self.equation)
        if eq.name != "f_mul":
            eq = TreeNode("f_mul", [eq])
        if eq.name == "f_mul":
            for child in eq.children:
                if child.name == "f_pow":
                    try:
                        n = int2(child.children[1].name[2:])
                        if n < 0:
                            for i in range(-n):
                                output.append(Eq("1")/Eq(child.children[0], "treeform"))
                        else:
                            for i in range(n):
                                output.append(Eq(child.children[0], "treeform"))
                    except:
                        output.append(Eq(child, "treeform"))
                else:
                    output.append(Eq(child, "treeform"))
        return output
    def __repr__(self):
        return string_equation(self.equation)

def substitute_val(eq, val):
    eq = replace(tree_form(eq.equation), tree_form("v_0"), tree_form("d_"+str(val)))
    return eq

def structure(eq, s, varlist={}):
    varlist = varlist.copy()
    def structure2(eq, s):
        nonlocal varlist
        if s.name[:2] == "u_":
            if s.name not in varlist.keys():
                varlist[s.name] = str_form(eq)
                return True
            elif varlist[s.name] == str_form(eq):
                return True
            return False
        if eq.name != s.name or len(eq.children)!=len(s.children):
            return False
        if eq.name in ["f_add", "f_mul"]:
            for item in itertools.permutations(list(s.children)):
                con = structure(eq.children[0], item[0], varlist.copy())
                if con is None:
                    continue
                fail = False
                for i in range(1,len(eq.children)):
                    con = structure(eq.children[i], item[i], con)
                    if con is None:
                        fail = True
                        break
                if not fail:
                    varlist = con
                    return True
            return False
        if any(not structure2(eq.children[i], s.children[i]) for i in range(len(eq.children))):
            return False
        return True
    
    if not structure2(eq, s):
        return None
    
    return varlist

def inversehandle(eq):
    if not hasattr(eq, "children") or not eq.children:
        return eq
    tmp = structure(copy.deepcopy(eq), tree_form('f_pow\n f_pow\n  u_0\n  u_1\n u_2'), {}.copy())
    if tmp is not None:
        for key in tmp.keys():
            tmp[key] = tree_form(tmp[key])
        exponent = TreeNode("f_mul", [tmp["u_1"], tmp["u_2"]])
        base = TreeNode("f_pow", [tmp["u_0"], exponent])
        return copy.deepcopy(base)
    if eq.name == "f_pow" and eq.children[0].name == "f_mul":
        exponent = copy.deepcopy(eq.children[1])
        ans = TreeNode("f_mul", [])
        for child in eq.children[0].children:
            ans.children.append(copy.deepcopy(TreeNode("f_pow", [child, exponent])))
        return ans
    return TreeNode(eq.name, [inversehandle(child) for child in eq.children])


def simplify(eq):
    if eq is None:
        return None
    if eq.name == "f_pow" and eq.children[0].name[:2] == "d_" and eq.children[1].name == "f_pow":
        
        if str_form(eq.children[1]) == "f_pow\n d_2\n d_-1":
            n = perfect_square_root(is_str_n(eq.children[0].name))
            if n is not None:
                return tree_form("d_"+str(n))
    if all(x.name[:2] == "d_" for x in eq.children):
        if eq.name == "f_sin":
            if eq.children[0].name == "d_0":
                return tree_form("d_0")
        if eq.name == "f_cos":
            if eq.children[0].name == "d_0":
                return tree_form("d_1")
        if eq.name == "f_pow":
            if int2(eq.children[0].name[2:]) == 0 and int2(eq.children[1].name[2:]) <= 0:
                return None
            if int2(eq.children[1].name[2:]) > 0:
                n = int2(eq.children[0].name[2:]) ** int2(eq.children[1].name[2:])
                return tree_form("d_"+str(n))
            if int2(eq.children[1].name[2:]) == 0:
                return None
            if abs(int2(eq.children[0].name[2:])) == 1 and int2(eq.children[1].name[2:]) == -1:
                return eq.children[0]
        elif eq.name == "f_add":
            arr = [int2(x.name[2:]) for x in eq.children]
            return tree_form("d_"+str(sum(arr)))
        elif eq.name == "f_mul":
            arr = [int2(x.name[2:]) for x in eq.children]
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
def replace_eq2(eq):
    
    if eq.name=="f_tan":
        return replace_eq2(TreeNode("f_div", [TreeNode("f_sin", [copy.deepcopy(eq.children[0])]), TreeNode("f_cos", [copy.deepcopy(eq.children[0])])]))

    arr = TreeNode(eq.name, [])
    for child in eq.children:
        arr.children.append(replace_eq2(child))
    return arr

def replace_eq(eq):
    if eq.name[2:] in ["tan", "sin", "cos"]:
        eq = TreeNode(eq.name, [tree_form(Eq(eq.children[0], "treeform").equation)])
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

def expand3(eq):
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

    return TreeNode(eq.name, [expand3(child) for child in eq.children])
def expand2(eq):
    
    if eq.name == "f_mul" or eq.name == "f_pow":
        if eq.name == "f_pow":
            eq = TreeNode("f_pow", [eq]) 
        ac = []
        addchild = []
        for child in eq.children:
            ac += [tree_form(x.equation) for x in Eq(child, "treeform").factor()]
        tmp3 = []
        for child in ac:
            tmp2 = []
            if child.name == "f_add":
                if child.children != []:
                    for child2 in child.children:
                        tmp2.append(Eq(child2, "treeform"))
                else:
                    tmp2 = [Eq(child, "treeform")]
            else:
                tmp3.append(Eq(child, "treeform"))
            if tmp2 != []:
                addchild.append(tmp2)
        tmp4 = Eq("1")
        for item in tmp3:
            tmp4 = tmp4 * item
        addchild.append([tmp4])
        def flatten(lst):
            flat_list = []
            for item in lst:
                if isinstance(item, list) and item == []:
                    continue
                if isinstance(item, list):
                    flat_list.extend(flatten(item))
                else:
                    flat_list.append(item)
            return flat_list
        
        if isinstance(addchild, list) and len(flatten(addchild))>0:
            add= Eq("0")
            for item in itertools.product(*addchild):
                mul = Eq("1")
                for item2 in item:
                    
                    mul = mul * item2
                add = add + mul
            eq = tree_form(add.equation)
    return TreeNode(eq.name, [expand2(child) for child in eq.children])
def numdem(equation):
    num = Eq("1")
    den = Eq("1")
    for item in equation.factor():
        t = tree_form(item.equation)
        if t.name == "f_pow" and int2(t.children[1].name[2:]) < 0:
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
            if len(new.children)==1:
                new = new.children[0]
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
    elif eq.name == "f_pow" and "v_" not in str_form(eq.children[1]):
        base, power = eq.children
        dbase = diff(Eq(base, "treeform"))
        b1 = Eq(power, "treeform") - Eq("1")
        bab1 = TreeNode("f_pow", [base, tree_form(b1.equation)])
        return Eq(power, "treeform")*Eq(bab1, "treeform")*dbase
def diffx2(equation, var="v_0"):
    if equation.name == "f_dif":
        if equation.children[0].name == var:
            return tree_form("d_1")
        return tree_form("d_0")
    return TreeNode(equation.name, [diffx2(child, var) for child in equation.children])
def diffany(equation):
    eq = tree_form(equation.equation)
    if eq.name == "f_dif":
        if eq.children[0].name[:2] == "v_":
            return Eq("1")
    return Eq(TreeNode(eq.name, [tree_form(diffany(Eq(child, "treeform")).equation) for child in eq.children]), "treeform")
def diffx(equation, var="v_0"):
    
    equation = diff(equation)
    
    equation = tree_form(equation.equation)
    equation = diffx2(equation, var)
    equation = Eq(equation, "treeform")
    
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
    try:
        equation = tree_form(solve(str_form(equation)))
    except:
        return None
        
    if equation is None:
        return None
    return Eq(equation, "treeform")
def lhospital(equation):
    
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
            g, h = diffx(n), diffx(d)
            equation = g/h
            return equation
        return None
    return None
def additionlimit(equation):
    equation = Eq(expand2(tree_form(equation.equation)), "treeform")

    final = Eq("0")
    equation = tree_form(equation.equation)
    if equation.name == "f_add":
        for child in equation.children:
            tmp = find_limit(Eq(child, "treeform"), 3, True)
            
            if tmp is not None:
                final += tmp
            else:
                final = None
                break
    else:
        final = None
    return final

def collect(eq):
    output = []
    def collect2(eq):
        if eq.name in ["f_pow", "f_sin", "f_cos"] and eq.children[0].name[:2] == "f_" and all(x not in str_form(eq.children[0]) for x in ["f_sin", "f_cos"]):
            
            output.append(str_form(eq.children[0]))
        for child in eq.children:
            collect2(child)
    collect2(eq)
    
    tmp = sorted(output, key=lambda x: len(x))[-1]
    tmp = solve(tmp)
    return tmp

def integratex(equation, var="v_0"):
    if diffx(equation, var) == Eq("0"):
        return Eq(var, "strform")*equation
    
    def varproc(term, var):
        if tree_form(diffx(term, var).equation).name[:2] == "d_":
            
            return diffx(term, var)
        return None
    const = [x for x in equation.factor() if var not in x.equation]
    tmp = Eq("1")
    for item in const:
        tmp = tmp * item
    const = tmp
    
    if const != Eq("1"):
        return const*integratex(equation/const, var)
    
    eq = tree_form(equation.equation)
    eq = expand2(eq)
    if eq.name == "f_add":
        return sum([integratex(Eq(child, "treeform"), var) for child in eq.children])
    if eq.name == var:
        return Eq(var, "strform")**Eq("2")/Eq("2")
    
    tmp4 = structure(eq, tree_form('f_pow\n u_0\n u_1'))
    if tmp4 is not None:
        tmp = varproc(Eq(tmp4["u_0"], "strform"), var)
        if tmp is not None:
            expo = Eq(tmp4["u_1"], "strform")+Eq("1")
            base = Eq(tmp4["u_0"], "strform")
            return base**expo / (expo*tmp)
    if eq.name == "f_cos":
        tmp = varproc(Eq(eq.children[0], "treeform"), var)
        if tmp is not None:
            return Eq(eq.children[0], "treeform").fx("sin")/tmp
    if eq.name == "f_sin":
        tmp = varproc(Eq(eq.children[0], "treeform"), var)
        if tmp is not None:
            return -1*Eq(eq.children[0], "treeform").fx("cos")/tmp
    if eq.name == "f_pow" and eq.children[0].name == "f_cos" and eq.children[1].name == "d_-2":
        tmp = varproc(Eq(eq.children[0].children[0], "treeform"), var)
        if tmp is not None:
            return Eq(eq.children[0].children[0], "treeform").fx("tan")/tmp
    if eq.name == "f_pow" and eq.children[0].name == "f_sin" and eq.children[1].name == "d_-2":
        tmp = varproc(Eq(eq.children[0].children[0], "treeform"), var)
        if tmp is not None:
            return 1/(Eq(eq.children[0].children[0], "treeform").fx("tan")*tmp)
        
    tmp4 = structure(eq, tree_form('f_mul\n f_pow\n  f_cos\n   u_0\n  d_-1\n f_sin\n  u_0'))
    
    if tmp4 is not None:
        tmp = varproc(Eq(tmp4["u_0"], "strform"), var)
        if tmp is not None:
            
            return -Eq(tmp4["u_0"], "strform").fx("cos").fx("abs").fx("log")/tmp

    tmp4 = structure(eq, tree_form('f_mul\n f_pow\n  f_sin\n   u_0\n  d_-1\n f_cos\n  u_0'))
    
    if tmp4 is not None:
        tmp = varproc(Eq(tmp4["u_0"], "strform"), var)
        if tmp is not None:
            return Eq(tmp4["u_0"], "strform").fx("sin").fx("abs").fx("log")/tmp

    tmp3= structure(eq, tree_form('f_mul\n f_sin\n  u_0\n f_pow\n  f_cos\n   u_0\n  d_-1\n f_sin\n  u_1\n f_pow\n  f_cos\n   u_1\n  d_-1\n f_sin\n  u_2\n f_pow\n  f_cos\n   u_2\n  d_-1'))
    if tmp3 is not None:
        for item in itertools.permutations(list(tmp3.keys())):
            item = [Eq(tmp3[x], "strform") for x in item]
            
            if item[0] + item[1] == item[2]:
                h = item[2].fx("tan")-item[1].fx("tan")-item[0].fx("tan")
                return integratex(Eq(replace_eq2(tree_form(h.equation)), "treeform"), var)
    if var == "v_0":
        x = Eq("y")-Eq(collect(eq), "strform")
        
        tmp3 = subs1(equation, x)
        if tmp3 is not None:
            return tmp3

def formula_1(eq):
    if eq.name == "f_pow" and eq.children[1].name[:2] == "d_" and abs(int(eq.children[1].name[2:])) == 2 and eq.children[0].name == "f_cos" and diffany2(Eq(eq.children[0].children[0], "treeform")).equation == "d_1":
        x = Eq(TreeNode("f_sin", [eq.children[0].children[0]]), "treeform")
        if int(eq.children[1].name[2:]) == 2:
            x = Eq("1")-x*x
        else:
            x = Eq("1")/(Eq("1")-x*x)
        return tree_form(x.equation)
    elif eq.name == "f_pow" and eq.children[1].name[:2] == "d_" and abs(int(eq.children[1].name[2:])) == 3 and eq.children[0].name == "f_cos" and diffany2(Eq(eq.children[0].children[0], "treeform")).equation == "d_1":
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
def find_limit(equation,depth=3,tail=False):
    
    ans = subslimit(equation)
    if ans is not None and "v_0" not in ans.equation:
        if "v_" in equation.equation.replace("v_0", ""):
            if "v_0" not in equation.equation:
                return equation
        else:
            return ans
    if depth == 0:
        return None
    
    for item in [lhospital,additionlimit,approx_limit]:
        ans = item(equation)
        
        if ans is not None:
            out = find_limit(ans, depth-1, tail)
            
            if out is not None:
                
                return out
    return None
def return_only_var(eq):
    output = []
    def helper(eq):
        nonlocal output
        if eq.name[:2] == "v_" and eq.name[2:].isdigit():
            output.append(eq.name)
        for child in eq.children:
            helper(child)
    helper(eq)
    output = list(set(output))
    if len(output)!=1:
        return None
    return output[0]

def const_mul(equation):
    output = Eq("d_1", "strform")
    
    for item in equation.factor():
        
        if "v_" not in item.equation:
            
            output = output * item
    return output

def degree(eq, var):
    
    if any("f_" in x in eq.equation for x in "sin cos tan".split(" ")):
        return 4
    count = 0
    while True:
        if eq == Eq("0") or count==4:
            break
        eq = diffx(eq, var)
        
        count += 1
    return count

def equal(a, b):
    x= Eq(a, "treeform")-Eq(b, "treeform")
    x = Eq(expand2(tree_form(x.equation)), "treeform")
    if x == Eq("0"):
        return True
    return False

def replace2(equation, find, r):
    if equal(equation, find):
        return r
    col = TreeNode(equation.name, [])
    for child in equation.children:
        col.children.append(replace2(child, find, r))
    return col

def contain(eq, term):
    if equal(eq, term):
        return True
    for child in eq.children:
        if contain(child, term):
            return True
    return False

def inverse(rhs,term):
    lhs = Eq("d_0", "strform")
    count = 6
    while not equal(rhs, term):
        if rhs.name == "f_add":
            for i in range(len(rhs.children)-1,-1,-1):
                if not contain(rhs.children[i], term):
                    lhs = lhs - Eq(rhs.children[i], "treeform")
                    rhs.children.pop(i)
        elif rhs.name == "f_mul":
            for i in range(len(rhs.children)-1,-1,-1):
                if not contain(rhs.children[i], term):
                    lhs = lhs / Eq(rhs.children[i], "treeform")
                    rhs.children.pop(i)
        elif rhs.name == "f_pow":
            lhs = lhs ** (Eq("1")/Eq(rhs.children[1], "treeform"))
            rhs = copy.deepcopy(rhs.children[0])
        if len(rhs.children) == 1:
            rhs = rhs.children[0]
        count -= 1
        if count == 0:
            return None
    return lhs


def inverse2(eq, var="v_0"):
    if isinstance(eq, Eq):
        eq = tree_form(eq.equation)
    tmp = inverse(eq, tree_form(var))
    return tmp


def subs1(equation, term):
    eq = tree_form(equation.equation)
    termeq = tree_form(term.equation)
    
    t = tree_form(inverse2(copy.deepcopy(termeq), "v_0").equation)
    t = expand2(t)
    eq = replace2(eq, tree_form("v_0"), t)
    
    eq2 = replace2(tree_form(diffx(inverse2(copy.deepcopy(termeq), "v_1")).equation), tree_form("v_0"), t)
    
    equation = Eq(eq, "treeform")/Eq(eq2, "treeform")
    
    equation = Eq(replace_eq(tree_form(equation.equation)), "treeform")
    equation = Eq(expand2(tree_form(equation.equation)), "treeform")

    tmp = integratex(equation, "v_1")
    
    if tmp is None:
        return None
    
    tmp = tree_form(tmp.equation)
    tmp = replace2(tmp, tree_form("v_1"), tree_form(inverse2(termeq, "v_1").equation))
    return Eq(tmp, "treeform")

def a2x2(eq):
    if eq.name == "f_pow" and str_form(eq.children[1]) == "f_pow\n d_-2\n d_-1" and eq.children[0].name == "f_add":
        equation = Eq(eq.children[0], "treeform")
        orig = equation
        for i in range(-2,3):
            equation = equation + Eq("d_"+str(i), "strform")
            eq2 = quadratic(tree_form(equation.equation))
            if eq2.name == "f_pow"and eq2.children[1].name == "d_2":
                var = Eq(eq2.children[0], "treeform")
                const = Eq("d_"+str(i), "strform")**Eq("1/2")
                t = (var-const)/(var+const)
                return t.fx("abs").fx("log") / (2*const)
            equation = orig
    return None
       
            
def quadratic(equation):
    output = None
    def quad(eq):
        
        tmp= return_only_var(eq)
        if tmp is not None:
            a,b,c = Eq("0"), Eq("0"), Eq("0")
            
            if eq.name != "f_add":
                eq = copy.deepcopy(TreeNode("f_add", [eq]))
            if eq.name == "f_add":
                for child in eq.children:
                    
                    cc = Eq(child, "treeform")
                    if degree(cc, tmp) == 3:
                        
                        a = const_mul(cc)
                    elif degree(cc, tmp) == 2:
                        
                        b = const_mul(cc)
                    elif degree(cc, tmp) == 1:
                        
                        c = const_mul(cc)
                    
            tmp = Eq(tmp, "strform")
            
            if a != Eq("0"):
                d = b**2 - 4*a*c
                if tree_form(d.equation).name[:2] == "d_":
                    n = int2(tree_form(d.equation).name[2:])
                    if n < 0:
                        return None
                else:
                    return None
                r1 = ( -b + (b**2 - 4*a*c)**Eq("1/2") )/(2*a)
                r2 = ( -b - (b**2 - 4*a*c)**Eq("1/2") )/(2*a)
                return (tmp-r1)*(tmp-r2)
            else:
                if b != Eq("0"):
                    return tmp+c/b
        return None

    if equation.name == "f_add":
        tmp = quad(copy.deepcopy(equation))
        if tmp is not None:
            return tree_form(tmp.equation)
    
    return TreeNode(equation.name, [quadratic(child) for child in equation.children])
def fraction2(equation):
    equation = tree_form(equation.equation)
    equation = formula_1(expand(common2(equation)))
    return Eq(equation, "treeform")
equation = None
while True:
    tmp = input(">>> ")
    try:
        orig = equation
        if tmp == "factor":   
            equation = Eq(quadratic(tree_form(equation.equation)), "treeform")
            
            print(equation)
        elif tmp == "expand":
             equation = Eq(expand2(tree_form(equation.equation)), "treeform")
             print(equation)
        elif tmp == "simplify":
             equation = Eq(simplify(tree_form(equation.equation)), "treeform")
             print(equation)
        elif tmp == "show":
            print(equation)
        elif tmp == "trig":
            equation = Eq(replace_eq(tree_form(equation.equation)), "treeform")
            print(equation)
        elif tmp == "fraction":
            equation = fraction2(equation)
            print(equation)
        elif tmp == "d/dx":
            equation = diffx(equation)
            print(equation)
        elif tmp[:2] == "Sd":
            tmp2 = integratex(equation, {"x":"v_0", "y":"v_1", "z":"v_2"}[tmp[-1]])
            if tmp2 is None:
                print("failed to integrate")
            else:
                equation = tmp2
                print(equation)
        elif tmp == "limit 0":
            
            print("limit x->0 " + str(equation) + " = " + str(find_limit(equation, 3, False)))
        else:
            equation = Eq(tmp,"stringequation")
            equation = Eq(replace_eq2(tree_form(equation.equation)), "treeform")
            print(equation)
    
    except Exception as error:
        equation = orig
        print(error)
