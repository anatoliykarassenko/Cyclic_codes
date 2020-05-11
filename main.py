import sympy
import itertools
from sympy.parsing.sympy_parser import parse_expr
from sympy.abc import x
import re
import numpy
from collections import Counter
from functools import reduce


def polinom(n):
    dl = list(sympy.factor(x ** n - 1, modulus=2).args)  # generating polynomial
    return dl


def Pol_generate(null_polinom, r):
    '''null_polinom outputs the null generating polynomial and all generating polynomials of degree r'''
    d2 = set()
    d3 = list()
    
    '''Loop for opening brackets'''
    for urav in null_polinom:
        if str(urav).find(')**') > 0:
            for i in range(int(str(urav)[-1])):
                new = parse_expr(str(urav)[1:-4])
                d3.append(new)
    if d3:
        null_polinom = d3
    for k_urav in range(1, len(null_polinom) + 1):
        for i in itertools.combinations(null_polinom, k_urav):  # iterating over bindings of polynomials
            z = reduce(lambda o, y: sympy.expand(o * y, modulus=2), i) # multiplication of polynomials
            d2.add(z)
    d2 = list(d2)
    
    if int(r) == 1:
        costill = list(d2[x] for x in range(0, int(r)))
        print(costill, 'All generating polynomials of degree r')
        return costill, costill
    else:
        ii = list(filter(lambda x: x[4:4 + len(r)] == r and x[4 + len(r)] == ' ',
                         str(d2).split(',')))  # search for a value with the required coefficient
        d_dopo = list(
            int(str(d2).split(',').index(ii[x])) for x in
            range(len(ii)))  # find the equation to the desired value
        
        
        D_por = list(d2[x] for x in d_dopo)
        print(D_por, 'All generating polynomials of degree r')
        
        return D_por[0], D_por


def matrixx(pp, n, r):
    '''Function for creating a generating matrix in two types: sympy and numpy. 
    Also finds the minimum code distance between matrix rows'''
    '''the first value is the generating polynomial, 
    the second value is the minimum code distance'''
    d3 = list()
    Pol_urav = list()
    xst = 1
    for i in range(n - int(r)):
        line_list = (re.sub(' ', '', str(sympy.expand(pp * xst))).split('+'))
        '''divide the equation by + and remove the spaces while shifting x'''
        # sorting from 1
        d3.append(sorted(line_list))
        Pol_urav.append(sympy.expand(pp * xst))
        xst *= x
    
    d_buf = list()
    
    for i in range(n - int(r)):
        x1st = 1
        now = []
        for ii in range(n):
            now.append(str(x1st))
            x1st *= x
        d_buf.append(now)
    
    for i in range(n - int(r)):
        d_plus = set(d_buf[i]) & set(d3[i])  # searching for elements and creating a matrix using sets
        # print(d_plus)
        ie = 0
        for ii in d_plus:
            d_buf[i][d_buf[i].index(ii)] = 1
        for ii in range(n):
            if d_buf[i][ie] != 1:
                d_buf[i][ie] = 0
                ie += 1
            else:
                ie += 1

    Pol_urav.append(0)  

    return numpy.array(d_buf)


def Research_max_dist(ishnpol, n, r):
    '''search function for a polynomial with the maximum code distance'''
    z1z = list(list(Spetr(ishnpol[x], n, r) for x in range(len(ishnpol))))
    d_mincod = list((list(z1z[x].keys())[1]) for x in range(len(ishnpol)))
    '''filling the array with the values of the minimum code distance for all 
    generating polynomials'''
    print(d_mincod, 'minimum code distance for polynomials of degree r')
    max_distance = ishnpol[d_mincod.index(max(d_mincod))]  
    return max_distance, z1z[d_mincod.index(max(d_mincod))]


def provepolinom(Porpolinom, n):
    provpol = sympy.factor((x ** n - 1) * Porpolinom ** -1, modulus=2)
    return sympy.expand(provpol, modulus=2)


def Spetr(Porpolinom, n, r):
    d_buf = list()
    x1st = 1
    # print(d3)
    for i in range(n - int(r)):
        d_buf.append(x1st)  # similarly to matrixx we collect a polynomial of degree r
        x1st *= x
    U = list()  # Empty list
    
    for i in range(1, n - int(r)):
        vrcomb = list(itertools.combinations(d_buf, i))  # search for all possible information vectors
        for ii in vrcomb:
            if str(ii)[-2::] == ',)':
                U.append(sympy.expand(parse_expr(str(ii)[1:-2]) * Porpolinom, modulus=2))  # with the number of elements 1
            else:
                zz1 = '+'.join(str(ii).split(','))  # a temp array
                U.append(sympy.expand(parse_expr(str(zz1)) * Porpolinom, modulus=2))
    U.append(sympy.expand(parse_expr('+'.join(str(d_buf)[1:-1].split(','))) * Porpolinom, modulus=2))
    '''Adding a product of a polynomial with information vectors with a number 
    of elements greater than 1 to the list'''
    xxx = list(str(str(ico).count('+') + 1) for ico in U)
    
    sspectr = {'0': 1}  # sorting the dictionary so that it starts with 0
    coount = Counter(xxx)
    sspectr.update(coount)  # unsorted spectrum
    ssss = dict()
    for key in sorted(sspectr, key=lambda i: int(i)):
        ssss.update({int(key): sspectr[key]})
    
    return ssss


print('Введите n')
nn = int(input())  # nn parameter
print('Введите r')
rr = input()
dl = polinom(nn)  # all polynomials with sympy
pp1, D_porvse = Pol_generate(dl, rr)
# print(pp1, D_porvse)
# print(por_mat[0])
maxcodpoly, spectral = Research_max_dist(D_porvse, nn, rr)  # search for a polynomial with the maximum code distance
print(maxcodpoly, 'generating polynomial')

por_mat = matrixx(maxcodpoly, nn, rr)  # creating a generating matrix
print(por_mat, ' generating matrix')  # output of the matrix with the maximum code distance
ProvP = provepolinom(maxcodpoly, nn)
print(ProvP, 'check polynomial')
print(spectral, 'Code spectrum')
