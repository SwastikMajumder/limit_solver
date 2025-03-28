import math
import copy
import itertools
from base import *
from fractions import Fraction
import parser
def int_nth_root(a, n):
    if a < 0 and n % 2 == 0:
        return None
    x = round(math.pow(a, 1 / n))
    if x**n == a:
        return x
    return None
def perfect_square_root(n, nn=2):
    root = None
    if n < 0:
        root = int_nth_root(-n,nn)
        if root is not None:
            root = -root
    else:
        root = int_nth_root(n,nn)
    if root is None:
        return None
    return root if root ** nn == abs(n) else None
def pf(n, nn=2):
    if not isinstance(n, Fraction):
        n = Fraction(n)
    a = n.numerator
    b = n.denominator
    if perfect_square_root(a, nn) is not None and perfect_square_root(b, nn) is not None:
        return Fraction(perfect_square_root(a, nn),perfect_square_root(b, nn))
    return None
def int2(string):
    tmp = Fraction(string)
    if tmp.denominator == 1:
        return tmp.numerator
    return tmp
def flatten_tree(node):
    if not node.children:
        return node
    if node.name in ("f_add", "f_mul", "f_and", "f_or"):
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
    
    while "f_sub" in str_form(eq) or "f_div" in str_form(eq):
        eq = a1(eq)
    eq = flatten_tree(eq)
    while "f_neg" in str_form(eq) or "f_inv" in str_form(eq):
        eq = a2(eq) 
    return eq
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
        
        if b.denominator != 1:
            tmp2 = pf(a, b.denominator)
            if tmp2 is None:
                return None
            a = tmp2
        try:
            a = a ** b.numerator
        except:
            return None
        return a
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
    if "f_list" in str_form(eq):
        return TreeNode(eq.name, [simple(copy.deepcopy(child)) for child in eq.children])
    if (eq.name[:2] == "f_" and eq.name != "f_add") or eq.name[:2] in ["d_","v_","s_"]:
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
            return TreeNode("f_mul", [expand2(a),c])
    arr = TreeNode(eq.name, [])
    for child in eq.children:
        arr.children.append(common(child))
    return arr
def structure(eq, s, varlist={}):
    varlist = varlist.copy()
    def structure2(eq, s):
        nonlocal varlist
        if s.name[:2] in ["u_", "p_"]:
            if "v_" in str_form(eq) and s.name[:2] == "p_":
                return False
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
def common2(eq):
    eq = common(eq)
    return TreeNode(eq.name, [common2(child) for child in eq.children])



def powermerge(eq):
    def base(eq, index):
        return eq.children[index].children[0] if eq.children[index].name == "f_pow" else eq.children[index]
    def expo(eq, index):
        return eq.children[index].children[1] if eq.children[index].name == "f_pow" else tree_form("d_1")
    if eq.name == "f_mul":
        change = True
        while change:
            change = False
            l = len(eq.children)
            for i in range(2, l + 1):
                for item in itertools.combinations(range(l), i):
                    if all(base(eq, item[0]) == base(eq, item2) for item2 in item):
                        merged_exponent = TreeNode("f_add", [expo(eq, item2) for item2 in item])
                        merged_term = base(eq, item[0]) ** merged_exponent
                        indices_to_remove = sorted(item, reverse=True)
                        for j in indices_to_remove:
                            eq.children.pop(j)
                        eq.children.append(merged_term)
                        if len(eq.children) == 1:
                            eq = copy.deepcopy(eq.children[0])
                        change = True
                        break
                if change:
                    break
    return TreeNode(eq.name, [powermerge(child) for child in eq.children])

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
    
    tmp = structure(eq, tree_form("s_e")**(tree_form("u_0")*tree_form("u_1").fx("log")))
    if tmp is not None:
        return tree_form(tmp["u_1"])**tree_form(tmp["u_0"])
    
    if eq.name == "f_log" and eq.children[0].name == "f_pow" and eq.children[0].children[0].name == "s_e":
        return inversehandle(eq.children[0].children[1])
    if eq.name == "f_pow" and eq.children[0].name == "s_e" and eq.children[1].name == "f_log":
        return inversehandle(eq.children[1].children[0])
    if (eq.name == "f_sin" and eq.children[0].name == "f_arcsin") or (eq.name == "f_cos" and eq.children[0].name == "f_arccos") or (eq.name == "f_tan" and eq.children[0].name == "f_arctan"):
        return inversehandle(eq.children[0].children[0])
    if (eq.name == "f_cos" and eq.children[0].name == "f_arcsin") or (eq.name == "f_sin" and eq.children[0].name == "f_arccos"):
        eq2 = eq.children[0].children[0]
        eq2 = (tree_form("d_1") - eq2*eq2)**(tree_form("d_1")/tree_form("d_2"))
        return inversehandle(eq2)
    if (eq.name == "f_arcsin" and eq.children[0].name == "f_sin") or (eq.name == "f_arccos" and eq.children[0].name == "f_cos") or (eq.name == "f_arctan" and eq.children[0].name == "f_tan"):
        return inversehandle(eq.children[0].children[0])
    
    return TreeNode(eq.name, [inversehandle(child) for child in eq.children])
def simplify(eq, lite=False):
    def iscont(a, b, c=None):
        if c is None:
            a, b = [tree_form(x).name for x in [a,b]]
        else:
            a, b, c = [tree_form(x).name for x in [a,b,c]]
        if (c is None and all(x[:2] == "d_" for x in [a,b])) or all(x[:2] == "d_" for x in [a,b,c]):
            
            if c is None:
                
                a, b = [int(x[2:]) for x in [a,b]]
                if a == 0:
                    return None
                a = Fraction(a,1)**Fraction(b,1)
                
                return frac2(tree_form("d_" + str(a)))
            else:
                a, b, c = [int(x[2:]) for x in [a,b,c]]
            
            s = "d_"
            
            if a < 0 and c == -1:
                return None
            if c == -1:
                a = pf(a, abs(b))
                if a is None:
                    return None
                if b < 0:
                    a = Fraction(1,1)/a
                return frac2(tree_form("d_" + str(a)))
        return None
    if eq is None:
        return None
    if not lite:
        tmp = structure(eq, tree_form('f_pow\n u_0\n f_pow\n  u_1\n  u_2'))
        if tmp is not None:
            
            tmp2 = iscont(*[tmp[x] for x in sorted(tmp.keys())])
            
            if tmp2 is not None:
                return tmp2
        
        tmp = structure(eq, tree_form('f_pow\n u_0\n u_1'))
        if tmp is not None:
            
            tmp2 = iscont(*[tmp[x] for x in sorted(tmp.keys())])
            
            if tmp2 is not None:
                return tmp2
    tmp = structure(eq, tree_form('f_mul\n d_1\n u_0'))
    if tmp is not None:
        return simplify(tree_form(tmp["u_0"]), lite)
    
    tmp = structure(eq, tree_form('f_add\n d_0\n u_0'))
    if tmp is not None:
        return simplify(tree_form(tmp["u_0"]), lite)

    tmp = structure(eq, tree_form('f_mul\n d_0\n u_0'))
    if tmp is not None:
        return tree_form("d_0")

    tmp = structure(eq, tree_form('f_pow\n u_0\n d_0'))
    if tmp is not None:
        return tree_form("d_1")
    
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
        child = simplify(child, lite)
        if child is None:
            return None
        arr.children.append(child)
    return arr

def dowhile(eq, fx):
    while True:
        orig = copy.deepcopy(eq)
        eq = copy.deepcopy(fx(eq))
        if eq is None:
            return None
        if eq == orig:
            return orig
def solve2(eq):
    eq = copy.deepcopy(eq)
    for item in [convert_sub2neg, expand_eq, flatten_tree, inversehandle, powermerge, simplify]:
        eq = dowhile(eq, item)
    return eq
def solve(eq):
    return dowhile(eq, solve2)

def factorgen(eq):
    output = []
    if eq.name != "f_mul":
        eq = TreeNode("f_mul", [eq])
    if eq.name == "f_mul":
        for child in eq.children:
            if child.name == "f_pow":
                if child.children[0].name == "s_e":
                    output.append(child)
                    continue
                try:
                    n = int2(child.children[1].name[2:])
                    if n < 0:
                        for i in range(-n):
                            output.append(tree_form("d_1")/child.children[0])
                    else:
                        for i in range(n):
                            output.append(child.children[0])
                except:
                    output.append(child)
            else:
                output.append(child)
    return [copy.deepcopy(simplifylite(x)) for x in output]
    
def substitute_val(eq, val, var="v_0"):
    eq = replace(eq, tree_form(var), tree_form("d_"+str(val)))
    return eq


def product_to_sum(eq):
    
    lst = factorgen(eq)
    
    if len(lst) == 1:
        return lst[0]

    if len(lst) == 2:
        a, b = lst
        
        if a.name == "f_sin" and b.name == "f_sin":
            return ((a.children[0] - b.children[0]).fx("cos") - (a.children[0] + b.children[0]).fx("cos")) / tree_form("d_2")
        elif a.name == "f_cos" and b.name == "f_cos":
            return ((a.children[0] - b.children[0]).fx("cos") + (a.children[0] + b.children[0]).fx("cos")) / tree_form("d_2")
        elif a.name == "f_sin" and b.name == "f_cos":
            return ((a.children[0] + b.children[0]).fx("sin") + (a.children[0] - b.children[0]).fx("sin")) / tree_form("d_2")
        elif a.name == "f_cos" and b.name == "f_sin":
            return ((a.children[0] + b.children[0]).fx("sin") - (a.children[0] - b.children[0]).fx("sin")) / tree_form("d_2")
        elif a.name == "f_tan" and b.name == "f_tan":
            return ((a.children[0].fx("sec") * b.children[0].fx("sec") - tree_form("d_1")) / 
                    (a.children[0].fx("sec") * b.children[0].fx("sec") + tree_form("d_1")))
        elif a.name == "f_sin" and b.name == "f_tan":
            return a.children[0].fx("sin") * b.children[0].fx("sec")
        elif a.name == "f_cos" and b.name == "f_tan":
            return a.children[0].fx("cos") * b.children[0].fx("sec")
    
    first, rest = lst[0], lst[1:]
    s = tree_form("d_0")
    eq = solve(expand2(first * product_to_sum(solve(TreeNode("f_mul", rest)))))
    if eq.name == "f_add":
        for child in eq.children:
            s += product_to_sum(child)
            s = solve(s)
    else:
        s = eq
    return s
def replace_eq3(eq):
    eq = product_to_sum(eq)
    eq = replace_eq2(eq)
    eq = solve(eq)
    return eq

def replace_eq2(eq):
    if eq.name=="f_tan":
        return replace_eq2(TreeNode("f_div", [TreeNode("f_sin", [copy.deepcopy(eq.children[0])]), TreeNode("f_cos", [copy.deepcopy(eq.children[0])])]))
    if eq.name == "f_sec":
        return replace_eq2(tree_form("d_1")/eq.children[0].fx("cos"))
    if eq.name == "f_cosec":
        return replace_eq2(tree_form("d_1")/eq.children[0].fx("sin"))
    if eq.name == "f_cot":
        return replace_eq2(eq.children[0].fx("cos")/eq.children[0].fx("sin"))


    def isneg(eq):
        eq = tree_form(eq)
        if eq.name[:2] != "d_":
            return False
        if int(eq.name[2:]) >= 0:
            return False
        return True
    
    tmp = structure(eq, (tree_form("u_0")*tree_form("p_0")).fx("sin"))
    if tmp is not None and isneg(tmp["p_0"]):
        return tree_form("d_-1")*(tree_form(tmp["u_0"])*tree_form(tmp["p_0"])*tree_form("d_-1")).fx("sin")

    tmp = structure(eq, (tree_form("u_0")*tree_form("p_0")).fx("cos"))
    if tmp is not None and isneg(tmp["p_0"]):
        return (tree_form(tmp["u_0"])*tree_form(tmp["p_0"])*tree_form("d_-1")).fx("cos")
    
    arr = TreeNode(eq.name, [])
    for child in eq.children:
        arr.children.append(replace_eq2(child))
    return arr

def replace_eq(eq):
    if eq.name[2:] in ["tan", "sin", "cos"]:
        eq = TreeNode(eq.name, [solve(eq.children[0])])
    
    tmp = structure(eq, (tree_form("u_0")*tree_form("d_-1")).fx("sin"))
    if tmp is not None:
        return tree_form("d_-1")*tree_form(tmp["u_0"]).fx("sin")
    
    if eq.name in ["f_sin", "f_cos"] and eq.children[0].name == "f_add" and len(eq.children[0].children)==2:
        if eq.name == "f_sin":  
            a = TreeNode("f_sin", [eq.children[0].children[0]])
            b = TreeNode("f_cos", [eq.children[0].children[1]])
            c = TreeNode("f_cos", [eq.children[0].children[0]])
            d = TreeNode("f_sin", [eq.children[0].children[1]])
            
            return replace_eq(a*b+c*d)
        elif eq.name == "f_cos":
            a = TreeNode("f_cos", [eq.children[0].children[0]])
            b = TreeNode("f_cos", [eq.children[0].children[1]])
            c = TreeNode("f_sin", [eq.children[0].children[0]])
            d = TreeNode("f_sin", [eq.children[0].children[1]])
            
            return replace_eq(a*b-c*d)

    tmp = structure(eq, tree_form('f_sin\n f_mul\n  d_2\n  u_0'))
    if tmp is not None:
        var = tree_form(tmp["u_0"])
        return solve(tree_form("d_2")*var.fx("sin")*var.fx("cos"))

    tmp = structure(eq, tree_form('f_cos\n f_mul\n  d_2\n  u_0'))
    if tmp is not None:
        var = tree_form(tmp["u_0"])
        return solve(var.fx("cos")*var.fx("cos")-var.fx("sin")*var.fx("sin"))
        
    arr = TreeNode(eq.name, [])
    for child in eq.children:
        arr.children.append(replace_eq(child))
    return arr

def expand(eq):
    if eq.name == "f_mul":
        addchild = [[child2 for child2 in child.children] for child in eq.children if child.name == "f_add"]
        otherchild = [child for child in eq.children if child.name != "f_add"]
        if len(otherchild) == 1 and otherchild[0].name[:2] == "d_" and len(addchild) == 1 and isinstance(addchild, list):
            add= tree_form("d_0")
            
            for item in addchild[0]:
                add += otherchild[0]*item
                
            return add
    return TreeNode(eq.name, [expand(child) for child in eq.children])

def addm(m1, m2):
    m1, m2 = copy.deepcopy(m1), copy.deepcopy(m2)
    if m1.name != "f_list":
        return copy.deepcopy(TreeNode("f_add", [m1, m2])), tree_form("d_0")
    coll = TreeNode("f_list", [])
    for i in range(len(m1.children)):
        
        coll.children.append(addm(m1.children[i], m2.children[i])[0])
    return coll, m2

def scale(m, number):
    if m.name != "f_list":
        return m*number
    return TreeNode("f_list", copy.deepcopy([scale(child, number) for child in m.children]))

def conv_mat(eq):
    if eq.name == "f_list":
        return [conv_mat(child) for child in eq.children]
    else:
        return eq

def conv_list(eq):
    if isinstance(eq, list):
        return TreeNode("f_list", [conv_list(child) for child in eq])
    else:
        return eq

def matrix_multiply(matrix1, matrix2):
    matrix1, matrix2= copy.deepcopy(matrix1), copy.deepcopy(matrix2)
    if len(matrix1[0]) != len(matrix2):
        return None

    output = []
    
    for i in range(len(matrix1)):
        element = []
        for j in range(len(matrix2[0])):
            coll = TreeNode("f_add", [])
            for k in range(len(matrix2)):
                coll.children.append(TreeNode("f_mul", [matrix1[i][k], matrix2[k][j]]))
            element.append(coll)
        output.append(element)
    return output

def transpose(matrix):
    return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

def convert_nested_list(lst):
    if isinstance(lst, int):
        return tree_form("d_"+str(lst))
    if isinstance(lst, list):
        return [convert_nested_list(item) for item in lst]
    return lst

def determinant(matrix):
    matrix = convert_nested_list(matrix)
    return determinant2(matrix)

def determinant2(matrix):
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    if n == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    det = tree_form("d_0")
    for col in range(n):
        sub_matrix = [row[:col] + row[col+1:] for row in matrix[1:]]
        det += tree_form("d_"+str(-1 ** col)) * matrix[0][col] * determinant2(sub_matrix)
    return det

def circumcenter(matrix):
    x1, y1 = matrix[0]
    x2, y2 = matrix[1]
    x3, y3 = matrix[2]
    a = x1**tree_form("d_2")+y1**tree_form("d_2")
    b = x2**tree_form("d_2")+y2**tree_form("d_2")
    c = x3**tree_form("d_2")+y3**tree_form("d_2")
    n1 = determinant([[a,y1,1],[b,y2,1],[c,y3,1]])
    n2 = determinant([[x1,a,1],[x2,b,1],[x3,c,1]])
    d = determinant([[x1,y1,1],[x2,y2,1],[x3,y3,1]])
    two = tree_form("d_2")
    return [n1/(two*d),n2/(two*d)]
def matrix_simp(eq):
    eq = copy.deepcopy(eq)
    if eq.name == "f_eq":
        return TreeNode(eq.name, [matrix_simp(eq.children[0]), tree_form("d_0")])
    if eq.name == "f_transpose":
        if eq.children[0].name != "f_list":
            eq.children[0] = matrix_simp(eq.children[0])
        eq.children[0] = conv_mat(eq.children[0])
        eq.children[0] = transpose(eq.children[0])
        return conv_list(eq.children[0])
    
    if eq.name == "f_circumcenter":
        if eq.children[0].name != "f_list":
            eq.children[0] = matrix_simp(eq.children[0])
        eq.children[0] = conv_mat(eq.children[0])
        eq.children[0] = circumcenter(eq.children[0])
        return conv_list(eq.children[0])
    
    if eq.name == "f_add":
        
        if eq.children[0].name != "f_list":
            eq.children[0] = matrix_simp(eq.children[0])
        collect = copy.deepcopy(eq.children[0])
        
        for i in range(1,len(eq.children)):
            if eq.children[i].name != "f_list":
                eq.children[i] = matrix_simp(eq.children[i])
            collect = addm(collect, eq.children[i])[0]
            
        return collect
    if eq.name == "f_mul" and any(child.name == "f_list" for child in eq.children):
        for i in range(len(eq.children)):
            if eq.children[i].name == "f_transpose":
                eq.children[i] = matrix_simp(eq.children[i])
        collect = copy.deepcopy([child for child in eq.children if child.name == "f_list"][0])
        collect = conv_mat(collect)
        index = [index for index, child in enumerate(eq.children) if child.name == "f_list"][0]
        
        for i in range(len(eq.children)):
            if i == index:
                continue
            if "f_list" in str_form(eq.children[i]):
                collect = matrix_multiply(collect, conv_mat(eq.children[i]))
            else:
                collect = scale(conv_list(collect), eq.children[i])
                collect = conv_mat(collect)
        
        return conv_list(collect)
    if eq.name == "f_pow":
        out = TreeNode("f_mul", [])
        for i in range(int(eq.children[1].name[2:])):
            out.children.append(copy.deepcopy(eq.children[0]))
        return matrix_simp(out)
    return eq
def simplifylite(eq):
    return simplify(eq, True)
def expand2(eq):
    eq = copy.deepcopy(eq)
    if eq.name == "f_mul" or eq.name == "f_pow":
        if eq.name == "f_pow":
            eq = copy.deepcopy(TreeNode("f_pow", [eq]) )
        ac = []
        addchild = []
        for child in eq.children:
            ac += factorgen(child)
        tmp3 = []
        for child in ac:
            tmp2 = []
            if child.name == "f_add":
                if child.children != []:
                    for child2 in child.children:
                        tmp2.append(child2)
                else:
                    tmp2 = [child]
            else:
                tmp3.append(child)
            if tmp2 != []:
                addchild.append(tmp2)
        tmp4 = tree_form("d_1")
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
            add= tree_form("d_0")
            for item in itertools.product(*addchild):
                mul = tree_form("d_1")
                for item2 in item:
                    mul = mul * item2
                    mul = dowhile(mul, simplifylite)
                add = add + mul
                add = dowhile(add, simplifylite)
            eq = copy.deepcopy(add)
            
    return TreeNode(eq.name, [expand2(child) for child in eq.children])

def numdem(equation):
    num = tree_form("d_1")
    den = tree_form("d_1")
    for item in factorgen(equation):
        t = item
        if t.name == "f_pow" and int2(t.children[1].name[2:]) < 0:
            den = den*item
        else:
            num = num*item
    return [solve(num), solve(tree_form("d_1")/den)]

def diff(eq):
    eq = copy.deepcopy(eq)
    eq = dowhile(eq, simplifylite)
    eq = solve(eq)
    if "v_" not in str_form(eq):
        return tree_form("d_0")
    if eq.name == "f_add":
        add = tree_form("d_0")
        for child in eq.children:
            add += diff(child)
        return add
    elif eq.name == "f_pow" and eq.children[0].name == "s_e":
        return diff(eq.children[1])*eq
    elif eq.name == "f_tan":
        return diff(eq.children[0])/(eq.children[0].fx("cos")*eq.children[0].fx("cos"))
    elif eq.name == "f_log":
        return diff(eq.children[0])*(tree_form("d_1")/eq.children[0])
    elif eq.name == "f_arcsin":
        return diff(eq.children[0])/(tree_form("d_1")-eq.children[0]*eq.children[0])**parser.take_input("1/2")
    elif eq.name == "f_arctan":
        return diff(eq.children[0])/(tree_form("d_1")+eq.children[0]*eq.children[0])
    elif eq.name == "f_mul":
        add = tree_form("d_0")
        for i in range(len(eq.children)):
            new = copy.deepcopy(eq)
            new.children.pop(i)
            if len(new.children)==1:
                new = new.children[0]
                
            add += diff(dowhile(eq.children[i], simplifylite))*new
        add = dowhile(add, simplifylite)
        return add
    elif eq.name == "f_sin":
        eq.name = "f_cos"
        return diff(eq.children[0])*eq
    elif eq.name == "f_cos":
        eq.name = "f_sin"
        return tree_form("d_-1")*diff(eq.children[0])*eq
    
    elif eq.name[:2] == "v_":
        return TreeNode("f_dif", [eq])
    
    elif eq.name == "f_pow" and "v_" not in str_form(eq.children[1]):
        base, power = eq.children
        dbase = diff(base)
        b1 = power - tree_form("d_1")
        
        bab1 = TreeNode("f_pow", [base, b1])
        
        return power * bab1 * dbase
    
def diffx2(equation, var="v_0"):
    if equation.name == "f_dif":
        if equation.children[0].name == var:
            return tree_form("d_1")
        return tree_form("d_0")
    return TreeNode(equation.name, [diffx2(child, var) for child in equation.children])

def diffany(eq):
    
    if eq.name == "f_dif":
        if eq.children[0].name[:2] == "v_":
            return tree_form("d_1")
    return TreeNode(eq.name, [diffany(child) for child in eq.children])

def diffx(equation, var="v_0"):
    equation = diff(equation)
    equation = diffx2(equation, var)
    return solve(equation)
def diffany2(equation):
    equation = diff(equation)
    equation = diffany(equation)
    equation = solve(equation)
    return equation
def approx(equation):
    con = True
    while con:
        con = False
        for item in factorgen(equation):
            if item.name == "f_sin":
                
                if diffx(diffx(item.children[0])) != tree_form("d_0"):
                    continue
                equation = equation/item
                equation = equation*item.children[0]
                con = True
                break
            elif item.name == "f_pow" and item.children[1].name == "d_-1" and item.children[0].name == "f_sin":
                if diffx(diffx(item.children[0].children[0])) != tree_form("d_0"):
                    continue
                equation = equation*item.children[0]
                equation = equation/item.children[0].children[0]
                con = True
                break
            elif item.name == "f_cos":
                if diffx(diffx(item.children[0])) != tree_form("d_0"):
                    continue
                equation = equation/item
                x = item.children[0]
                equation = equation*(tree_form("d_1")-x*x/tree_form("d_2"))
                con = True
                break
            elif item.name == "f_pow" and item.children[1].name == "d_-1" and item.children[0].name == "f_cos":
                if diffx(diffx(item.children[0].children[0])) != tree_form("d_0"):
                    continue
                equation = equation*item.children[0]
                x = item.children[0].children[0]
                equation = equation/(tree_form("d_1")-x*x/tree_form("d_2"))
                con = True
                break
    equation = solve(equation)
    
    return equation

def subslimit(equation):
    equation = expand2(equation)
    equation = simplify(equation)
    if equation is None:
        return None
    
    equation = substitute_val(equation, 0)
    equation = dowhile(equation, simplify)
    try:
        equation = solve(equation)
    except:
        return None
    if equation is None:
        return None
    return equation
def lhospital(equation):
    equation = simplify(equation)
    if equation is None:
        return None
    
    e = substitute_val(equation, 0)
    e = dowhile(e, simplify)
    if e is None:
        n, d = numdem(equation)
        ans1 = subslimit(n)
        ans2 = subslimit(d)
        if ans1 is not None and ans2 is not None and ans1 == tree_form("d_0") and ans2 == tree_form("d_0"):
            g, h = diffx(n), diffx(d)
            equation = g/h
            return solve(equation)
        return None
    return None
def additionlimit(equation):
    equation = expand2(equation)
    final = tree_form("d_0")
    
    if equation.name == "f_add":
        for child in equation.children:
            tmp = find_limit(child, 3, True)
            if tmp is not None:
                final += tmp
                final = solve(tmp)
            else:
                final = None
                break
    else:
        final = None
    return final

def collec(eq, var):
    output = []
    def collect2(eq):
        if eq.name in ["f_pow", "f_sin", "f_cos"] and eq.children[0].name[:2] != "v_" and var in str_form(eq.children[0]):
            output.append(str_form(eq.children[0]))
        if eq.name == "f_pow" and eq.children[0].name == "s_e" and "v_" in str_form(eq):
            if eq.children[1].name[:2] != "v_":
                output.append(str_form(eq.children[1]))
            output.append(str_form(eq))
        for child in eq.children:
            collect2(child)
    collect2(eq)
    tmp = sorted(output, key=lambda x: len(x))
    tmp = [solve(tree_form(x)) for x in tmp]
    return tmp
constant_term= 0
def ccc():
    global constant_term
    constant_term += 1
    return tree_form("v_c"+str(constant_term-1))
def integratex(equation, var="v_0"):
    
    if var not in ["v_0", "v_1", "v_2"]:
        return None
    
    if not contain(equation, tree_form(var)):
        return tree_form(var)*equation
    
    def varproc(term, var):
        if diffx(term, var).name[:2] == "d_":
            return diffx(term, var)
        return None
    eq2 = copy.deepcopy(equation)
    
    const = [x for x in factorgen(eq2) if var not in str_form(x)]
    
    tmp = tree_form("d_1")
    for item in const:
        tmp = tmp * item
    const = tmp
    const = solve(const)
    
    if const != tree_form("d_1"):
        r = integratex(solve(equation/const), var)
        if r is not None:
            return const*r
    
    eq = solve(equation)

    eq2 = solve(expand2(eq))
    
    if eq.name in ["f_mul", "f_pow"] and eq2.name == "f_add":
        tmp4 = integratex(eq2, var)
        if tmp4 is not None:
            return tmp4

    eq3 = solve(expand2(replace_eq3(eq)))
    
    if solve(expand2(eq3-eq2)) != tree_form("d_0") and solve(expand2(eq3-eq)) != tree_form("d_0"):
        tmp4 = integratex(eq3, var)
        if tmp4 is not None:
            return tmp4
        
    if eq.name == "f_add":
        success = True
        s = integratex(eq.children[0], var)
        if s is not None:
            for child in eq.children[1:]:
                u = integratex(child, var)
                if u is None:
                    success = False
                    break
                s += u
            if success:
                return s
    
    if eq.name == var:
        return tree_form(var)**tree_form("d_2")/tree_form("d_2")
    
    if eq.name == "f_pow" and eq.children[0].name == "s_e":
        
        tmp = varproc(eq.children[1], var)
        
        if tmp is not None:
            
            return eq/tmp
    if eq.name == "f_log" and eq.children[0].name == var:
        return eq*tree_form(var)-tree_form(var)
    
    tmp4 = structure(eq, tree_form('f_pow\n u_0\n u_1'))
    if tmp4 is not None:
        tmp = varproc(tree_form(tmp4["u_0"]), var)
        if tmp is not None:
            
            expo = solve(tree_form(tmp4["u_1"])+tree_form("d_1"))
            base = tree_form(tmp4["u_0"])
            
            if expo == tree_form("d_0"):
                return base.fx("abs").fx("log")/tmp
            return base**expo / (expo*tmp)
        
    if eq.name == "f_cos":
        tmp = varproc(eq.children[0], var)
        if tmp is not None:
            return eq.children[0].fx("sin")/tmp
    if eq.name == "f_sin":
        tmp = varproc(eq.children[0], var)
        if tmp is not None:
            return tree_form("d_-1")*eq.children[0].fx("cos")/tmp
        
    if eq.name == "f_pow" and eq.children[0].name == "f_cos" and eq.children[1].name == "d_-2":
        tmp = varproc(eq.children[0].children[0], var)
        if tmp is not None:
            return eq.children[0].children[0].fx("tan")/tmp
        
    if eq.name == "f_pow" and eq.children[0].name == "f_sin" and eq.children[1].name == "d_-2":
        tmp = varproc(eq.children[0].children[0], var)
        if tmp is not None:
            return tree_form("d_-1")/(eq.children[0].children[0].fx("tan")*tmp)

    tmp4 = structure(eq, tree_form('f_mul\n f_sin\n  u_0\n f_pow\n  f_cos\n   u_0\n  d_-2'))
    if tmp4 is not None:
        tmp = varproc(tree_form(tmp4["u_0"]), var)
        if tmp is not None:
            return tree_form("d_1")/(tree_form(tmp4["u_0"]).fx("cos")*tmp)

    tmp4 = structure(eq, tree_form('f_mul\n f_cos\n  u_0\n f_pow\n  f_sin\n   u_0\n  d_-2'))
    if tmp4 is not None:
        tmp = varproc(tree_form(tmp4["u_0"]), var)
        if tmp is not None:
            return tree_form("d_-1")/(tree_form(tmp4["u_0"]).fx("sin")*tmp)
    
    tmp4 = structure(eq, tree_form('f_mul\n f_pow\n  f_cos\n   u_0\n  d_-1\n f_sin\n  u_0'))
    if tmp4 is not None:
        tmp = varproc(tree_form(tmp4["u_0"]), var)
        if tmp is not None:
            return tree_form("d_-1")*tree_form(tmp4["u_0"]).fx("cos").fx("abs").fx("log")/tmp
    
    tmp4 = structure(eq, tree_form('f_mul\n f_pow\n  f_sin\n   u_0\n  d_-1\n f_cos\n  u_0'))
    if tmp4 is not None:
        tmp = varproc(tree_form(tmp4["u_0"]), var)
        if tmp is not None:
            return tree_form(tmp4["u_0"]).fx("sin").fx("abs").fx("log")/tmp

    tmp3 = structure(eq, tree_form('f_mul\n u_0\n f_pow\n  f_add\n   u_1\n   f_pow\n    u_2\n    u_3\n  f_mul\n   u_4\n   f_pow\n    u_5\n    d_-1'))
    if tmp3 is not None:
        tmp = [tree_form(tmp3[x]) for x in sorted(tmp3.keys(), key=lambda x: x)]
        if tmp[5] == tmp[3]:
            return integratex(tmp[2]**tmp[4]*tmp[0]*(tmp[1]+tmp[2]**(tree_form("d_-1")*tmp[3]))**(tmp[4]/tmp[3]), var)
    
    v2 = "v_"+str(int(var[2:])+1)
    for item in collec(eq, var):
        
        x = tree_form(v2)-item
        tmp3 = subs1(equation, x, var, v2)
        
        if tmp3 is not None:
            return tmp3
        
    return None
def formula_1(eq):
    if eq.name == "f_pow" and eq.children[1].name[:2] == "d_" and abs(int(eq.children[1].name[2:])) == 2 and eq.children[0].name == "f_cos" and diffany2(eq.children[0].children[0]).name == "d_1":
        x = TreeNode("f_sin", [eq.children[0].children[0]])
        if int(eq.children[1].name[2:]) == 2:
            x = tree_form("d_1")-x*x
        else:
            x = tree_form("d_1")/(tree_form("d_1")-x*x)
        return x
    
    elif eq.name == "f_pow" and eq.children[1].name[:2] == "d_" and abs(int(eq.children[1].name[2:])) == 3 and eq.children[0].name == "f_cos" and diffany2(eq.children[0].children[0]).name == "d_1":
        x = TreeNode("f_sin", [eq.children[0].children[0]])
        y = TreeNode("f_cos", [eq.children[0].children[0]])
        if int(eq.children[1].name[2:]) == 3:
            x = y-x*x*y
        else:
            x = tree_form("d_1")/(y-x*x*y)
        return x
    return TreeNode(eq.name, [formula_1(child) for child in eq.children])
def approx_limit(equation):
    orig = equation
    equation = approx(equation)
    if orig == equation:
        return None
    return equation
def find_limit(equation,depth=3,tail=False):
    ans = subslimit(equation)
    if ans is not None and "v_0" not in str_form(ans):
        if "v_" in str_form(equation).replace("v_0", ""):
            if "v_0" not in str_form(equation):
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
    output = tree_form("d_1")
    for item in factorgen(equation):
        if "v_" not in str_form(item):
            output = output * item
    return output
def degree(eq, var):
    if any("f_" in x in str_form(x) for x in "sin cos tan".split(" ")):
        return 4
    count = 0
    while True:
        if eq == tree_form("d_0") or count==4:
            break
        eq = diffx(eq, var)
        count += 1
    return count
def equal(a, b):
    x = a-b
    x = expand2(x)
    x = solve(x)
    if x == tree_form("d_0"):
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
def clearv(eq):
    
    if eq.name == "f_eq":
        equation=eq.children[0]
        arr = [x for x in factorgen(equation) if "v_" in str_form(x)]
        for i in range(len(arr)-1,-1,-1):
            if arr[i].name == "f_pow" and arr[i].children[1].name[:2] == "d_" and int(arr[i].children[1].name[2:]) < 0:
                arr.pop(i)
        m = arr[0]
        for item in arr[1:]:
            m = m * item
        return TreeNode(eq.name, [m, tree_form("d_0")])
    return TreeNode(eq.name, [clearv(child) for child in eq.children])
def fraction2(equation):
    equation = formula_1(expand(common2(equation)))
    equation = solve(equation)
    return equation
def inverse(rhs,term):
    lhs = tree_form("d_0")
    count = 6
    while not equal(rhs, term):
        
        if rhs.name == "f_add":
            
            if all(term in factorgen(child) for child in rhs.children):
                
                newrhs = solve(expand2(rhs/term))
                if not contain(newrhs, term):
                    rhs = term * newrhs
            else:
                for i in range(len(rhs.children)-1,-1,-1):
                    if not contain(rhs.children[i], term):
                        lhs = lhs - rhs.children[i]
                        rhs.children.pop(i)
        elif rhs.name == "f_mul":
            for i in range(len(rhs.children)-1,-1,-1):
                if not contain(rhs.children[i], term):
                    lhs = lhs / rhs.children[i]
                    rhs.children.pop(i)
        elif rhs.name == "f_pow" and contain(rhs.children[0], term):
            lhs = lhs ** (tree_form("d_1")/rhs.children[1])
            rhs = copy.deepcopy(rhs.children[0])
        elif rhs.name == "f_sin" and contain(rhs.children[0], term):
            lhs = lhs.fx("arcsin")
            rhs = copy.deepcopy(rhs.children[0])
        elif rhs.name == "f_arcsin" and contain(rhs.children[0], term):
            lhs = lhs.fx("sin")
            rhs = copy.deepcopy(rhs.children[0])
        elif rhs.name == "f_arccos" and contain(rhs.children[0], term):
            lhs = lhs.fx("cos")
            rhs = copy.deepcopy(rhs.children[0])
        elif rhs.name == "f_cos" and contain(rhs.children[0], term):
            lhs = lhs.fx("arccos")
            rhs = copy.deepcopy(rhs.children[0])
        elif rhs.name == "f_log" and contain(rhs.children[0], term):
            lhs = tree_form("s_e")**lhs
            rhs = copy.deepcopy(rhs.children[0])
        elif rhs.name == "f_pow" and rhs.children[0].name == "s_e" and contain(rhs.children[1], term):
            lhs = lhs.fx("log")
            rhs = copy.deepcopy(rhs.children[1].fx("log"))
        elif rhs.name == "f_tan" and contain(rhs.children[0], term):
            lhs = lhs.fx("arctan")
            rhs = copy.deepcopy(rhs.children[0])
        elif rhs.name == "f_arctan" and contain(rhs.children[0], term):
            lhs = lhs.fx("tan")
            rhs = copy.deepcopy(rhs.children[0])
        if len(rhs.children) == 1:
            rhs = rhs.children[0]
        count -= 1
        if count == 0:
            return None
        
    return solve(lhs)

def inverse2(eq, var="v_0"):
    
    def newfx(eq):
        eq = expand(common2(eq))
        eq = clearv(TreeNode("f_eq", [eq, tree_form("d_0")])).children[0]
        eq = solve(eq)
        return eq
    eq = dowhile(eq, newfx)
    
    eq2 = inverse(eq, tree_form(var))
    
    if eq2 is None:
        return None
    tmp = solve(expand2(eq2))
    return tmp

def subs1(equation, term, v1, v2):
    
    equation = solve(equation)
    eq = equation
    termeq = term
    
    t = inverse2(copy.deepcopy(termeq), v1)
    g = inverse2(termeq, v2)

    if g is None:
        return None
    
    if t is None:
        eq2 = diffx(g, v1)

        
        equation = solve(eq/eq2)
        
        equation = replace2(equation, g, tree_form(v2))
        
        if v1 in str_form(equation):
            return None
    else:
        t = expand2(t)
        
        eq = replace2(eq, tree_form(v1), t)
        eq2 = replace2(diffx(g, v1), tree_form(v1), t)
        
        
        equation = eq/eq2
        
        equation = solve(equation)
        
    lst = [ replace_eq(equation), equation]
    if solve(lst[0]-lst[1]) == tree_form("d_0"):
        lst = [lst[0]]
    for eq in lst:
        if v1 in str_form(eq):
            continue
        eq = solve(eq)
        
        tmp = integratex(eq, v2)
        
        if tmp is None:
            continue

        tmp = replace2(tmp, tree_form(v2), g)
        return tmp
    
    return None

def rref(matrix):
    rows, cols = len(matrix), len(matrix[0])
    lead = 0
    for r in range(rows):
        if lead >= cols:
            return matrix
        i = r
        while solve(matrix[i][lead]) == tree_form("d_0"):
            i += 1
            if i == rows:
                i = r
                lead += 1
                if lead == cols:
                    return matrix
        matrix[i], matrix[r] = matrix[r], matrix[i]
        lv = matrix[r][lead]
        matrix[r] = [solve(m / lv) for m in matrix[r]]
        for i in range(rows):
            if i != r:
                lv = matrix[i][lead]
                matrix[i] = [solve(m - lv * n) for m, n in zip(matrix[i], matrix[r])]
        lead += 1
    return matrix

def islinear(eq):
    
    eq =solve(eq)
    if eq.name == "f_pow" and "v_" in str_form(eq):
        return False
    for child in eq.children:
        out = islinear(child)
        if not out:
            return out
    return True

def linear(eqlist):
    final = []
    extra = []
    
    for i in range(len(eqlist)-1,-1,-1):
        if eqlist[i].name == "f_mul" and not islinear(expand2(eqlist[i])):
            if "v_" in str_form(eqlist[i]):
                eqlist[i] = TreeNode("f_mul", [child for child in eqlist[i].children if "v_" in str_form(child)])
                
            if all(islinear(child) for child in eqlist[i].children):
                for child in eqlist[i].children:
                    extra.append(TreeNode("f_eq", [child, tree_form("d_0")]))
                eqlist.pop(i)
            else:
                final.append(TreeNode("f_eq", [eqlist[i], tree_form("d_0")]))
                eqlist.pop(i)
    if extra != []:
        final.append(TreeNode("f_or", extra))
    if eqlist == []:
        if len(final)==1:
            return final[0]
        return TreeNode("f_and", final)
    eqlist = [eq for eq in eqlist if "v_" in str_form(eq)]
    
    if not all(islinear(eq) for eq in eqlist):
        
        return TreeNode("f_and", copy.deepcopy(final+eqlist))
    vl = []
    def varlist(eq):
        nonlocal vl
        if eq.name[:2] == "v_":
            vl.append(eq.name)
        for child in eq.children:
            varlist(child)
    
    for eq in eqlist:
        varlist(eq)
    vl = list(set(vl))
    if len(vl) > len(eqlist):
        return TreeNode("f_and", final+eqlist)
    m = []
    for eq in eqlist:
        s = copy.deepcopy(eq)
        
        row = []
        for v in vl:
            row.append(diffx(eq, v))
            s = replace(s, tree_form(v), tree_form("d_0"))
        row.append(s)
        m.append(row)
    
    m = rref(m)
    output = []
    for index, row in enumerate(m):
        count = 0
        for item in row[:-1]:
            if item == tree_form("d_1"):
                count += 1
                if count == 2:
                    break
            elif item == tree_form("d_0") and count == 1:
                break
        if count == 0:
            continue
        output.append(tree_form(vl[index])+row[-1])
    return TreeNode("f_and", final+[TreeNode("f_eq", [x, tree_form("d_0")]) for x in output])


def quadratic(equation):
    output = None
    def quad(eq):
        tmp= return_only_var(eq)
        if tmp is not None:
            eq = expand2(eq)
            eq = solve(eq)
            if eq.name == "f_add":
                for i in range(len(eq.children)):
                    if eq.children[i].name[:2] != "d_" and eq.children[i].name != "f_mul":
                        eq.children[i] = tree_form("d_1")*eq.children[i]
                if all(child.name[:2] != "d_" for child in eq.children):
                    eq = eq + tree_form("d_0")
                
                tmp = structure(eq, tree_form('f_add\n f_mul\n  f_pow\n   u_0\n   d_2\n  p_0\n f_mul\n  u_0\n  p_1\n p_2'))
                if tmp is not None:
                    a, b, c = [tree_form(x) for x in [tmp["p_0"], tmp["p_1"], tmp["p_2"]]]
                    var = tree_form(tmp["u_0"])
                    r1 = ( tree_form("d_-1")*b + (b**tree_form("d_2") - tree_form("d_4")*a*c)** (tree_form("d_1")/tree_form("d_2") ) )/(tree_form("d_2")*a)
                    r2 = ( tree_form("d_-1")*b - (b**tree_form("d_2") - tree_form("d_4")*a*c)** (tree_form("d_1")/tree_form("d_2") ) )/(tree_form("d_2")*a)
                    r1, r2 = solve(r1), solve(r2)
                    return a*solve((var-r1)*(var-r2))
                tmp = structure(eq, tree_form('f_add\n f_mul\n  f_pow\n   u_0\n   d_2\n  p_0\n p_1'))
                if tmp is not None:
                    a, b, c = [tree_form(x) for x in [tmp["p_0"], "d_0", tmp["p_1"]]]
                    var = tree_form(tmp["u_0"])
                    r1 = ( tree_form("d_-1")*b + (b**tree_form("d_2") - tree_form("d_4")*a*c)** (tree_form("d_1")/tree_form("d_2") ) )/(tree_form("d_2")*a)
                    r2 = ( tree_form("d_-1")*b - (b**tree_form("d_2") - tree_form("d_4")*a*c)** (tree_form("d_1")/tree_form("d_2") ) )/(tree_form("d_2")*a)
                    r1, r2 = solve(r1), solve(r2)
                    return a*solve((var-r1)*(var-r2))
                tmp = structure(eq, tree_form('f_add\n f_mul\n  u_0\n  p_0\n p_1'))
                if tmp is not None:
                    a, b, c = [tree_form(x) for x in ["d_0", tmp["p_0"], tmp["p_1"]]]
                    var = tree_form(tmp["u_0"])
                    return b*solve(var+c/b)
        return None
    if equation.name == "f_add":
        tmp = quad(copy.deepcopy(equation))
        if tmp is not None:
            if not isinstance(tmp, TreeNode):
                tmp = tmp[0] * tmp[1]
            return tmp
    return TreeNode(equation.name, [quadratic(child) for child in equation.children])

def eqhandle(eq):
    if eq.name == "f_eq":
        if "f_list" in str_form(eq):
            if eq.children[1] != tree_form("d_0"):
                eq.children[0] = eqhandle(eq.children[1] - eq.children[0])
                eq.children[1] = tree_form("d_0")
            else:
                eq.children[0] = eqhandle(eq.children[0])
            return eq
        return TreeNode(eq.name, [eqhandle(eq.children[1]-eq.children[0]), tree_form("d_0")])
    return TreeNode(eq.name, [eqhandle(copy.deepcopy(child)) for child in eq.children])
def rmeq(eq):
    if eq.name == "f_eq":
        return rmeq(eq.children[0])
    return TreeNode(eq.name, [rmeq(child) for child in eq.children])
def mat0(eq):
    def findeq(eq):
        out = []
        if eq.name != "f_list":
            return [str_form(eq)]
        else:
            for child in eq.children:
                out += findeq(child)
        return out
    eqlist = findeq(eq)
    eqlist = [tree_form(x) for x in eqlist]
    eqlist = [rmeq(x) for x in eqlist]
    eqlist = [TreeNode("f_mul", factorgen(x)) for x in eqlist if x != tree_form("d_0")]
    eqlist = [x.children[0] if len(x.children) == 1 else x for x in eqlist]
    return linear(copy.deepcopy(eqlist))
def and0(eq):
    if eq.name == "f_and":
        eq2 = copy.deepcopy(eq)
        eq2.name = "f_list"
        return mat0(eq2)
    elif eq.name == "f_eq":
        return mat0(eq)
    return TreeNode(eq.name, [and0(child) for child in eq.children])
def mat00(eq):
    if eq.name == "f_eq":
        if "f_list" in str_form(eq):
            return mat0(eq.children[0])
        return TreeNode(eq.name, [eq.children[0], tree_form("d_0")])
    return TreeNode(eq.name, [mat00(child) for child in eq.children])
def normal(string, variable=None):
    
    equation = None
    if variable is not None:
        equation = parser.take_input(string,list(variable.keys()))
    else:
        equation = parser.take_input(string,None)

    def varlist(eq):
        out = []
        if eq.name[:2] == "v_":
            out.append(eq.name)
        for child in eq.children:
            out += varlist(child)
        return sorted(list(set(out)), key=lambda x: int(x[2:]))
    
    def req(eq, key, variable, vl):
        if eq.name == "f_"+key:
            tmp = variable[key]
            for i in range(len(eq.children)):
                tmp = replace(copy.deepcopy(tmp), tree_form(vl[i]), eq.children[i])
            return tmp
        return TreeNode(eq.name, [req(child, key, variable, vl) for child in eq.children])

    
    if variable is not None:
        for key in variable.keys():
            equation = req(equation, key, variable, varlist(variable[key]))
    
    equation = eqhandle(equation)
    
    #equation = replace_eq2(equation)
    
    if "f_list" in str_form(equation):
        equation = copy.deepcopy(convert_sub2neg(equation))
        
        equation = matrix_simp(equation)
       
        equation = matrix_simp(equation) 

    return equation

def formula_2(eq, var="v_0"):
    orig = copy.deepcopy(eq)
    eq = replace(eq, tree_form(var).fx("sin")**tree_form("d_2"), (tree_form("d_1") - tree_form(var).fx("cos")**tree_form("d_2")))
    eq = flatten_tree(eq)
    eq = solve(expand2(eq))
    if "f_sin" not in str_form(eq):
        return eq
    eq = orig
    orig = copy.deepcopy(eq)
    eq = replace(eq, tree_form(var).fx("cos")**tree_form("d_2"), (tree_form("d_1") - tree_form(var).fx("sin")**tree_form("d_2")))
    eq = flatten_tree(eq)
    eq = solve(expand2(eq))
    if "f_cos" not in str_form(eq):
        return eq
    return eq

intconst = ["v_3", "v_4", "v_5", "v_6", "v_7", "v_8"]
def intreg(eq):
    global intconst
    if eq.name == "f_integrate":
        v = tree_form(intconst.pop(0))
        return integratex(eq.children[0], eq.children[1].name) + v
    return TreeNode(eq.name, [intreg(child) for child in eq.children])

equation = None
variable = None
while True:
    tmp = input(">>> ")
    if True:
        orig = equation
        if tmp == "factor":   
            equation = quadratic(equation)
            equation = clearv(equation)
            equation = simplify(equation)
            print(equation)
        elif tmp == "expand":
             equation = expand2(equation)
             equation = solve(equation)
             print(equation)
        elif tmp == "simplify":
             equation = clearv(equation)
             equation = solve(equation)
             print(equation)
        elif tmp == "show":
            print(equation)
        elif tmp == "trig0":
            equation = replace_eq2(equation)
            print(equation)
        elif tmp == "trig":
            equation = replace_eq(equation)
            print(equation)
        elif tmp == "trig1":
            equation = replace_eq3(equation)
            print(equation)
        elif tmp == "square":
            if equation.name == "f_eq":
                
                if equation.children[0].name == "f_add":
                    out = []
                    for i in range(len(equation.children[0].children)-1,-1,-1):
                        if equation.children[0].children[i].name == "f_pow" and str_form(equation.children[0].children[i].children[1]) == 'f_pow\n d_2\n d_-1':
                            out.append(equation.children[0].children[i])
                            equation.children[0].children.children.pop(i)
                            break
                    for item in out:
                        equation.children[1] += item
                    
                equation = TreeNode(equation.name, [expand2(equation.children[0]*equation.children[0]), expand(equation.children[1]*equation.children[1])])
            equation = solve(equation)
            equation = eqhandle(equation)
            print(equation)
        elif tmp == "fraction":
            
            equation = fraction2(equation)
            
            print(equation)
        elif tmp == "d/dx":
            equation = diffx(equation)
            print(equation)
        elif tmp[:2] == "Sd":
            tmp2 = integratex(equation, {"x":"v_0", "y":"v_1", "z":"v_2", "c": "v_5"}[tmp[-1]])

            if tmp2 is None:
                print("failed to integrate")
            else:
                equation = tmp2
                print(equation)
        elif tmp.split(" ")[0] == "func":
            if variable is None:
                variable = {}
            out = normal(copy.deepcopy(tmp.split(" ")[2]), None)
            out = intreg(out)
            variable[tmp.split(" ")[1]] = out
            print(out)
        elif tmp == "mat":
            if "f_list" in str_form(equation):
                
                equation = mat00(equation)
            else:
                equation = and0(equation)
            print(equation)
        elif tmp.split(" ")[0] == "limit":
            n = int(tmp.split(" ")[1])
            eq2 = TreeNode("f_add",[tree_form("v_0"), tree_form("d_"+str(n))])
            eq = replace2(equation, tree_form("v_0"), eq2)
            orig = equation
            equation = find_limit(eq, 3, False)
            print(f"limit x->{n} " + str(orig) + " = " + str(equation))
        else:
            out = normal(tmp, copy.deepcopy(variable))
            out = intreg(out)
            out = replace(out, tree_form("v_11"), tree_form("d_-1")**(tree_form("d_2")**tree_form("d_-1")))
            equation = solve(out)
            print(equation)
    '''     
    except Exception as error:
        equation = orig
        print(error)
    '''
