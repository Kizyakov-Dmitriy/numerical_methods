# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 10:42:02 2021

@author: Дима
"""

from numpy import linalg as LA
from numpy import inf
from random import randint, random
from time import time
import csv
from fractions import Fraction
import math
import copy
import pandas as pd

def trans(mat):
    return [[mat[i][j] for i in range(len(mat))]for j in range(len(mat))]

def add(mat1, mat2):
    return [[mat1[i][j] + mat2[i][j] for j in range(len(mat1))]for i in range(len(mat1))]

def sub(mat1, mat2):
    return [[mat1[i][j] - mat2[i][j] for j in range(len(mat1))]for i in range(len(mat1))]

def mult(mat1, mat2):
    ln = len(mat1) 
    return [[sum([mat1[i][k] * mat2[k][j] for k in range(ln)]) for j in range(ln)] for i in range(ln)]

def const(mat, c):
    c = float(c) if "j" not in c else complex(c) 
    return [[mat[i][j] * c for j in range(len(mat))]for i in range(len(mat))]

def alt(dic):
    rez = ""
    for i in str(len(dic)):
        rez += chr(int(i) + 1040)
    return rez

def lst_conv(srt):
    srt = srt.replace(" ", "")
    srt = srt.replace("i", "j")
    for i in range(len(srt)):
        if srt[i] == "j":
            l = srt.rfind("(", 0, i)
            r = srt.find(")", i)
            srt = srt[:l] + "[" + srt[l+1:r] +"]" + srt[r+1:]
    br_count = 0
    n = -1
    while len(srt) - n > 1:
        n += 1
        if srt[n] == "[":
            br_count += 1
            continue
        if srt[n] == "]":
            br_count -= 1
            continue
        if not br_count and srt[n].isdigit():
            k = n
            srt += " "
            while srt[k].isdigit() or srt[k] == ".":
                k += 1
            srt = srt[:n] + "[" + srt[n:k] + "]" + srt[k:]
            br_count += 1
            srt = srt[:-1]
            
    br_count = 0
    n = -1
    while len(srt) - n > 1:

        n += 1
        if srt[n] == "{":
            br_count += 1
            continue
        if srt[n] == "}":
            br_count -= 1
            continue
        if not br_count and srt[n].isupper():
            k = n
            srt += " "
            while srt[k].isupper():
                k += 1
            srt = srt[:n] + "{" + srt[n:k] + "}" + srt[k:]
            br_count += 1
            srt = srt[:-1]
    lst = []
    br_count = 0
    for i in range(len(srt)):
        if srt[i] in "}]":
            br_count -= 1
            continue
        if srt[i] in "{[":
            lst.append(srt[i+1:srt.find(chr(ord(srt[i])+2), i)])
            br_count += 1
        if not br_count:
            lst.append(srt[i])
    return lst

def lst_counter(lst, dic):
    while ")" in lst:
        for i, el in enumerate(lst):
            if el == "(":
                left_ = i
            if el == ")":
                right_ = i
                break
        a = lst_counter(lst[left_+1:right_], dic)
        lst = lst[:left_] + a[0] + lst[right_+1:]
        dic = a[1]
    
    while "t" in lst:
        for i, el in enumerate(lst):
            if el == "t":
                name = f"X{alt(dic)}"
                dic[name] = trans(dic[lst[i+1]])
                lst[i+1] = name
        for _ in range(lst.count("t")):
            lst.remove("t")
    while "*" in lst:
        for i, el in enumerate(lst):
            if el == "*":

                name = f"X{alt(dic)}"
                if lst[i-1][0].isdigit():
                    dic[name] = const(dic[lst[i+1]], lst[i-1])
                    lst[i] = name
                    lst[i-1] = ""
                    lst[i+1] = ""
                    break
                elif lst[i+1][0].isdigit():
                    dic[name] = const(dic[lst[i-1]], lst[i+1])
                    lst[i] = name
                    lst[i-1] = ""
                    lst[i+1] = ""
                    break
                else:
                    dic[name] = mult(dic[lst[i-1]], dic[lst[i+1]])
                    lst[i] = name
                    lst[i-1] = ""
                    lst[i+1] = ""
                    break
        for _ in range(lst.count("")):
            lst.remove("")
    while "+" in lst or "-" in lst:
        for i, el in enumerate(lst):
            if el == "+":
                name = f"X{alt(dic)}"
                dic[name] = add(dic[lst[i+1]], dic[lst[i-1]])
                lst[i] = name
                lst[i-1] = ""
                lst[i+1] = ""
                break
            if el == "-":
                name = f"X{alt(dic)}"
                dic[name] = add(dic[lst[i-1]], dic[lst[i+1]])
                lst[i] = name
                lst[i-1] = ""
                lst[i+1] = ""
                break
        for _ in range(lst.count("")):
            lst.remove("")
    
    return [lst, dic]

def calc(srt, dic):
    """Функция вычисляющая матричные выражения по строке и словарю
            Параметры:
                    srt (str): Строка-выражение содержащая сложение(+), вычитание(-), умножение(*), транспонирование(t(mat))
                    dic (dict): словарь {название матрицы : матрица}
            Возвращаемое значение: 
                    конечная матрица"""
    a = lst_counter(lst_conv(srt), dic)
    return a[1][a[0][0]]

def recurs(mat):
    if len(mat[0]) == 1:
        return mat[0][0]
    k = 0
    for i in range(1, len(mat)):
        if mat[i].count(0) > mat[k].count(0):
            k = i
    det = 0
    for ind, el in enumerate(mat[k]):
        det += (-1)**(k+ind) * el * recurs([[mat[i][j] for j in range(len(mat)) if j != ind] for i in range(len(mat)) if i != k])
    return det

def det_check(mat):
    norm_inf = max([sum([abs(i) for i in row]) for row in mat])
    try:
        inv = LA.inv(mat)
    except LA.LinAlgError:
        assert 1==2, "Матрица вырожденная"
    norm_inf_inv = max([sum([abs(i) for i in row]) for row in inv])
    assert norm_inf * norm_inf_inv < 500, "Матрица слабо обусловленна"
    assert len(mat) == len(mat[0]), "Матрица не является правильной"
    return mat

def determ(mat):
    """Функция вычисления определителя матрицы с проверкой на вырожденность и слабую обусловленность
            Параметры: 
                    mat (list): Матрица в виде списка списков
            Возвращаемое значение: определитель матрицы либо исклюсительная ситуация"""
    return recurs(det_check(mat))

def mat_print(x): 
    print(pd.DataFrame(data=x))

def obusl(mat, inv):
    norm_inf = max([sum([abs(i) for i in row]) for row in mat])
    norm_inf_inv = max([sum([abs(i) for i in row]) for row in inv])
    return norm_inf * norm_inf_inv

def Yakobi(mat, toch = 0.001):
    mat_inp = copy.deepcopy(mat)
    X0 = [0] * len(mat)
    for i in range(len(mat)):
        zn = mat[i][i]
        mat[i][i] = 0
        for j in range(len(mat) + 1):
            if j < len(mat):
                mat[i][j] = (-1) * mat[i][j] / zn
            else:
                mat[i][j] = mat[i][j] / zn
    c = 0
    X1 = []
    while True:
        c = c + 1
        for i in range(len(mat)):
            X1.append(sum([mat[i][j] * X0[j] for j in range(len(mat))]) + mat[i][len(mat)])
        if max([abs(X1[i] - X0[i]) for i in range(len(X0))]) < toch:
            break
        if c > 1000:
            X1 = "матрица расходится"
            mat_inp = [[1, 2, 0], [2, 3.9999, 0]]
            break
        else:
            X0 = X1
            X1 = []
    
    mat_rev = LA.inv(list(map(lambda x: x[:-1], mat_inp))).tolist()
    
    return {"mat" : mat_inp, "inv" : mat_rev, "sol" : X1}

def Jordan_dig(mat):
    ln = len(mat)
    mat_inp = copy.deepcopy(mat)
    for i in range(ln):
        mat[i] += [1 if i == j else 0 for j in range(ln)]
    for i in range(ln):
        for j in range(ln):
            for k in range(2*ln, i-1, -1):
                mat[i][k] /= mat[i][i]
            if i != j:
                for k in range(2 * ln, i-1, -1):
                    mat[j][k] -= mat[i][k] * mat[j][i]
    mat_rev = []
    rez = []
    for i in range(ln):
        mat_rev.append(mat[i][ln+1:])
        rez.append(mat[i][ln])
    return {"mat" : mat_inp, "inv" : mat_rev, "sol" : rez} 

fr = lambda x: Fraction(str(x))

def to_float(mat):
    rez = []
    for i in range(len(mat)):
        rez.append(list(map(lambda x: float(x[:x.find("/")]) / float(x[x.find("/") + 1:]), mat[i])))
        return rez

def Jordan_fr(mat, ratio = True):
    ln = len(mat)
    mat_inp = copy.deepcopy(mat)
    
    for i in range(ln):
        mat[i] += [1 if i == j else 0 for j in range(ln)]
    
    for i in range(ln):
        for j in range(2*ln):
            mat[i][j] = fr(mat[i][j])
    
        
    for i in range(ln):
        for j in range(ln):
            for k in range(2*ln, i-1, -1):
                mat[i][k] /= mat[i][i]
            if i != j:
                for k in range(2 * ln, i-1, -1):
                    mat[j][k] -= mat[i][k] * mat[j][i]
    mat_rev = []
    rez = []
    for i in range(ln):
        if ratio:
            mat_rev.append(list(map(lambda x: f"{x.numerator}/{x.denominator}", mat[i][ln+1:])))
            rez.append(f"{mat[i][ln].numerator}/{mat[i][ln].denominator}")
        else:
            mat_rev.append(list(map(float, mat[i][ln+1:])))
            rez.append(float(mat[i][ln]))
    return {"mat" : mat_inp, "inv" : mat_rev, "sol" : rez}

def SLAU(matrix):
    '''
    Последовательно вычисляет решение СЛАУ различными методами и предоставляет решение и обратную матрицу основываясь на устойчивости решения
            Параметры:
                    matrix (list): матрица в виде списка списков размера n|n+1, где последний из столбцов - столбец свободных членов
            Возвращаемое значение:
                    none
    '''
    comp = complex == type(matrix[0][0])
    Ya = Yakobi(copy.deepcopy(matrix))
    Ya_ob = obusl(list(map(lambda x: x[:-1], Ya["mat"])), Ya["inv"])
    if Ya_ob < 100:
        print("Матрица была посчитана методом Якоби \n")
        print("Введенная матрица")
        mat_print(Ya["mat"])
        print("Обратная матрица")
        mat_print(Ya["inv"])
        print("Вектор решений матрицы")
        mat_print(Ya["sol"])
        print("Число обусловленности полученной матрицы:", Ya_ob)
        return None
    Jo_dig = Jordan_dig(copy.deepcopy(matrix))
    Jo_dig_ob = obusl(list(map(lambda x: x[:-1], Jo_dig["mat"])), Jo_dig["inv"])
    if Jo_dig_ob < 100 or comp:
        print("Матрица была посчитана методом Жордана-Гаусса \n")
        print("Введенная матрица")
        mat_print(Jo_dig["mat"])
        print("Обратная матрица")
        mat_print(Jo_dig["inv"])
        print("Вектор решений матрицы")
        mat_print(Jo_dig["sol"])
        print("Число обусловленности полученной матрицы:", Jo_dig_ob)
        print("Число обусловленности полученной матрицы при использовании метода Якоби:", Ya_ob)
        return None
    if not comp:
        Jo_fr = Jordan_fr(copy.deepcopy(matrix))
        Jo_fr_ob = obusl(list(map(lambda x: x[:-1], Jo_fr["mat"])), to_float(Jo_fr["inv"]))
        print("Матрица была посчитана методом Жордана-Гаусса в правильных дробях\n")
        print("Введенная матрица")
        mat_print(Jo_fr["mat"])
        print("Обратная матрица")
        mat_print(Jo_fr["inv"])
        print("Вектор решений матрицы")
        mat_print(Jo_fr["sol"])
        print("Число обусловленности полученной матрицы:", Jo_fr_ob)
        print("Число обусловленности полученной матрицы при использовании метода Якоби:", Ya_ob)
        print("Число обусловленности полученной матрицы при использовании метода Жардана-Гаусса:", Jo_dig_ob)
        return None
 
def CalculateEnergy(sequence, matrix):

    n = len(sequence)
    E = 0
    for i in range(n-1):
        E += Metric(sequence[i], sequence[i+1], matrix)
        
    
    E += Metric(sequence[-1], sequence[0], matrix)
    
    return E

def Metric(A,B, matrix):
    distance = matrix[A][B]
    return distance

def GenerateStateCandidate(seq):
    n = len(seq)               #определяем размер последовательности
    i = randint(0,n-1)      #генерируем целое случайное число
    j = randint(0,n-1)

    

    if i > j:
        seq = seq[:j] + seq[i:j:-1] + [seq[j]] + seq[i+1:]#обращаем подпоследовательность
    else:
        seq = seq[:i] + seq[j:i:-1] + [seq[i]] + seq[j+1:]#обращаем подпоследовательность
        
    return seq

def preobraz(seq, indexes):
    
    rezult = []
    for i in indexes:
        rezult.append([seq[i]])
    
    return rezult

def GetTransitionProbability(E, T):
    P = math.exp(-E/T)
    return P

def MakeTransit(probability):
    if (probability > 1) or (probability < 0):
        raise ValueError('Violation of argument constraint')
        
    value = random()

    if value <= probability:
        a = True
    else:
        a = False 
        
    return a

def DecreaseTemperature(initialTemperature, k):
    T = initialTemperature * 0.1 / k
    return T

def SimulatedAnnealing(n, matrix):
    '''
    Метод оптимизации алгоритмом отжига
            Параметры:
                    n (int): количество городов
                    matrix (list): весовая матрица в формате списка списков чисел типа float 
                    
            Возвращаемое значение:
                    (list(list, float)): первый элемент списка последовательность обхода городов, второй - общий путь
    '''
    initialTemperature = 10
    endTemperature = 0.00001
    state = [i for i in range(n)]   #начальное состояние
    
    currentEnergy = CalculateEnergy(state, matrix) #вычисляем энергию для первого состояния
    T = initialTemperature
    
    k = 0 
    for i in range(1,100000):  #на всякий случай ограничеваем количество итераций     

        stateCandidate = GenerateStateCandidate(state) #получаем состояние-кандидат - в индексах
        #Переводим наши индексы обратно в последовательность точек
#         stateCandidate = preobraz(state, stateCandidate)
        
        candidateEnergy = CalculateEnergy(stateCandidate, matrix)  #вычисляем его энергию
        
        #Считаем кол-во повторов
        if candidateEnergy == currentEnergy:
            k += 1
        elif candidateEnergy < currentEnergy:
            k = 0
        
        if candidateEnergy < currentEnergy:  #если кандидат обладает меньшей энергией
            currentEnergy = candidateEnergy  #то оно становится текущим состоянием
            state = stateCandidate
        else:
            p = GetTransitionProbability(candidateEnergy-currentEnergy, T) # иначе, считаем вероятность
            if MakeTransit(p): #и смотрим, осуществится ли переход
                currentEnergy = candidateEnergy
                state = stateCandidate

        T = DecreaseTemperature(initialTemperature, i) #уменьшаем температуру
        
        if (T <= endTemperature) or (k == 100): #условие выхода
            #Вывод:
            stroka_rez = ''
            for el in state:
                stroka_rez += f'{el+1} -> '
            stroka_rez += str(state[0]+1)
            #print('Текущий вариант: ', stroka_rez, 'Стоимость: ', currentEnergy)
            print(f'Выход произошел на {k}ой одинаковой иттерации после {i} иттераций')
            return [[el + 1 for el in state] + [state[0]+1], currentEnergy]
            
            break
        #Вывод:
        stroka_rez = ''
        for el in state:
            stroka_rez += f'{el+1} -> '
        stroka_rez += str(state[0]+1)
    else:
        return [[el + 1 for el in state] + [state[0]+1], currentEnergy]
        #print('Текущий вариант: ', stroka_rez, 'Стоимость: ', currentEnergy)

def ACO(N, mat, a = 0.6, b = 0.65, ro = 0.7, it = 1000):
    '''
    Метод оптимизации алгоритмом отжига
            Параметры:
                    n (int): количество городов
                    matrix (list): весовая матрица в формате списка списков чисел типа float 
                    a (float): константа отвечаюая за жадность алгоритме (меньше - более жадный)
                    b (float): иная костанта альтернатива a
                    ro (float): часть феромона, остающаяся после каждой итерации
                    it (float): количество итераций максимально совершаемых алгоритмом

            Возвращаемое значение:
                    (list(list, float)): первый элемент списка последовательность обхода городов, второй - общий путь
    '''
    mat = [[mat[i][j] if mat[i][j] != "inf" else inf for j in range(N)]for i in range(N)] # Заменяем "inf" на бесконечность из numpy
    tau = 100 #Начальный уровень феромона
    Q = 100
    fer = [[tau for i in range(N)]for j in range(N)] #Матрица феромона
    way = []
    length = inf
    
    q = 0
    out = False
    
    for i in range(it):
        fer1 = [[fer[i][j] * ro for j in range(N)]for i in range(N)]
        
        for k in range(N):
            X = [k]
            now = X[-1]
            
            while len(X) != N:
                
                ver = []
                for ind, t in enumerate(mat[now]):
                    if ind not in X:
                        ver.append(fer[now][ind]**a * (1/mat[now][ind])**b)
                ver = [q / sum(ver) for q in ver]
                ch = random()
                j = 0
                sm = 0
                while sm < ch:
                    j += 1
                    sm += ver[j-1]
                k = -1
                for ind, t in enumerate(mat[now]):
                    if ind not in X:
                        k += 1
                    if k == j:
                        break
                
                X.append(k)
                now = k
                
                
            c = 0
            X.append(X[0])
            for i1 in range(N):
                c += mat[X[i1]][X[i1+1]]

            if c < length:
                length = c
                way = X
                q = 0
            else:
                q += 1

            
            for i1 in range(N):
                fer1[X[i1]][X[i1+1]] += Q/c  
            if q == it*N*0.2:
                print(f'Выход произошел на {q}ой одинаковой иттерации после {i} иттераций')
                out = True
                break
        fer = fer1.copy()
        if out:
            break
    way = list(map(lambda x: x+1, way))
    return [way, length]