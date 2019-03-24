import math
import sympy
from sympy import sqrt, Mul, Add, Pow, simplify


class Symbol:

    @staticmethod
    def err_add(v1, e1, v2, e2):
        return sqrt(e1 ** 2 + e2 ** 2)

    @staticmethod
    def err_mul(v1, e1, v2, e2):
        return v1 * v2 * sqrt((e1 / v1) ** 2 + (e2 / v2) ** 2)

    @staticmethod
    def err_pow(v1, e1, v2):
        return (v1 ** v2) * (e1 / v1) * abs(v2)

    @staticmethod
    def is_number(n):
        t = type(n)
        if t == int or t == float:
            return True
        elif isinstance(n, Symbol):
            return False
        else:
            # sympy object in this case
            try:
                return n.is_number
            except:
                return False

    @staticmethod
    def from_mean(name, *args):
        mean, err = Symbol.std_err_of_mean(*args)
        return Symbol(name, val=mean, err=err)

    @staticmethod
    def std_err_of_mean(*args):
        n = len(args)
        # args = [sympy.Number(arg) for arg in args]
        mean = sum(args) / len(args)
        temp = 0
        for arg in args:
            temp += (arg - mean) ** 2
        return mean, math.sqrt(temp / (n - 1)) / math.sqrt(n)

    @staticmethod
    def sqrt(v):
        return Symbol(sym=sqrt(v.sym),
                      err_sym=sympy.Rational(1, 2) * v.err_sym,
                      val_set=v.val_set.copy())

    def __init__(self, sym="", err_sym="", val=None, err=None, val_set=dict()):

        self.val_set = val_set

        if type(sym) == str:
            self.sym = sympy.Symbol(sym)

            if self.is_number(val):
                self.val_set[self.sym] = val
            else:
                self.val_set[self.sym] = self.sym
        else:
            self.sym = sym

        if type(err_sym) == str:
            if len(err_sym) == 0:
                self.err_sym = sympy.Symbol('sigma_' + sym)
            else:
                self.err_sym = sympy.Symbol(err_sym)

            if self.is_number(err):
                self.val_set[self.err_sym] = err
            else:
                self.val_set[self.err_sym] = self.err_sym
        else:
            self.err_sym = err_sym

    def __add__(self, other):
        val_set = self.val_set.copy()

        if self.is_number(other):
            return Symbol(sym=self.sym + other,
                          err_sym=self.err_sym,
                          val_set=val_set)
        else:
            val_set.update(other.val_set)
            return Symbol(sym=self.sym + other.sym,
                          err_sym=self.err_add(
                              self.sym, self.err_sym, other.sym, other.err_sym),
                          val_set=val_set)

    def __radd__(self, other):
        val_set = self.val_set.copy()

        if self.is_number(other):
            return Symbol(sym=self.sym + other,
                          err_sym=self.err_sym,
                          val_set=val_set)

        else:
            raise NotImplementedError(
                str(other) + " is not a supported oprand")

    def __mul__(self, other):
        val_set = self.val_set.copy()

        if self.is_number(other):
            return Symbol(sym=self.sym * other,
                          err_sym=self.err_sym * other,
                          val_set=val_set)
        else:
            val_set.update(other.val_set)
            return Symbol(sym=self.sym * other.sym,
                          err_sym=self.err_mul(
                              self.sym, self.err_sym, other.sym, other.err_sym),
                          val_set=val_set)

    def __rmul__(self, other):
        val_set = self.val_set.copy()

        if self.is_number(other):
            return Symbol(sym=self.sym * other,
                          err_sym=self.err_sym * other,
                          val_set=val_set)
        else:
            raise NotImplementedError(
                str(other) + " is not a supported oprand")

    def __pow__(self, other):
        val_set = self.val_set.copy()

        if self.is_number(other):
            if other == 0.5:
                return self.sqrt(self)
            else:
                return Symbol(sym=self.sym ** other,
                              err_sym=self.err_pow(
                                  self.sym, self.err_sym, other),
                              val_set=val_set)
        else:
            raise NotImplementedError(
                str(other) + " is not a supported oprand")
            # print('Warning: raise to a power that has uncertainty')
            # val_set.update(other.val_set)
            # return Symbol(sym=self.sym ** other.sym,
            #               err_sym=self.err_pow(
            #                   self.sym, self.err_sym, other.sym),
            #               val_set=val_set)

    def __neg__(self):
        return Symbol(sym=-self.sym, err_sym=self.err_sym, val_set=self.val_set.copy())

    def __truediv__(self, other):
        return self * (other ** -1)

    def __rtruediv__(self, other):
        # x / Symbol = x * Symbol^-1
        return other * (self ** -1)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        # x - Symbol = (-Symbol) + x
        return (-self) + other

    def __str__(self):
        val, err = self.eval()
        return "{} | {}\n{} | {}".format(val, err, str(self.sym), str(self.err_sym))

    def simp(self):
        self.sym = sympy.simplify(self.sym)
        self.err_sym = sympy.simplify(self.err_sym)

    def eval(self):
        symbols = list(self.val_set.items())
        val = sympy.N(self.sym.subs(symbols))
        err = sympy.N(self.err_sym.subs(symbols))
        return val, err

    def eval_sym(self):
        # convert number to string to disable evaluation
        symbols_no_eval = [(k, str(v)) for (k, v) in self.val_set.items()]
        val_no_eval = self.sym.subs(symbols_no_eval)
        err_no_eval = self.err_sym.subs(symbols_no_eval)
        return val_no_eval, err_no_eval

    def expr(self):
        return self.sym, self.err_sym

    def display(self, name=""):
        from IPython.display import display, Math, Markdown

        val, err = self.eval()
        val_sym, err_sym = self.eval_sym()
        options = {
            'long_frac_ratio': 10
        }

        content = """$$
        \\begin{{align*}}
            \\textbf{{Numeric: }} \\\\
                              {0} &= \\ {1} \\pm {2} \\\\
            \\textbf{{Symbolic: }}& \\\\
            \\text{{Value: }} {0} &= {3} \\\\
                                  &= {4} \\\\
            \\text{{Error: }} \\sigma_{{{0}}}&= {5} \\\\
                                  &= {6} \\\\
        \\end{{align*}}
            $$""".format(name, val, err, sympy.latex(self.sym), sympy.latex(val_sym), sympy.latex(self.err_sym),  sympy.latex(err_sym))

        display(Math(content))


if __name__ == "__main__":

    a = Symbol('a', val=10, err=0.1)
    b = Symbol('b', val=100, err=2)
    c = Symbol('c', val=150, err=3)
    d = Symbol('d', val=123, err=1)

    f = 1 + a + 2 * b * c - d / 2
    print(f)
