# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 13:31:44 2021

@author: Дима
"""
from numpy import arange, linalg, abs
import math
import copy
import matplotlib.pyplot as plt
import pywt
import sys
from scipy.fft import rfft, rfftfreq


def fact(n):
    if n == 1 or n == 0:
        return n
    k = 1
    for i in range(2, n + 1):
        k = k * i
    return k


def Lagrandg(lst):
    """Интерполирование методом Лагранжа
            Параметры:
                    lst (list): список точек в формате [x i ,y i ]
            Возвращаемое значение: 
                    rez (list): список точек в формате [x i , y i , f i ] (fi – значение интерполированной функции в точках)
                    """
    rez = []
    for i in range(1, len(lst)):
        for x in arange(lst[i-1][0], lst[i][0], (lst[i][0] - lst[i-1][0])/10):
            summ = 0
            for j in range(len(lst)):
                mult = lst[j][1]
                for n in range(len(lst)):
                    if j != n:
                        mult *= (x - lst[n][0]) / (lst[j][0] - lst[n][0])
                summ += mult
            
            rez.append([x, summ])
    summ = 0
    for j in range(len(lst)):
        mult = lst[j][1]
        for n in range(len(lst)):
            if j != n:
                mult *= (lst[-1][0] - lst[n][0]) / (lst[j][0] - lst[n][0])
        summ += mult
    rez.append([lst[-1][0], summ])
    
    
    
    plt.plot([i[0] for i in lst], [i[1] for i in lst], "o")
    plt.plot([i[0] for i in rez], [i[1] for i in rez], "r")
    plt.title('Интерполяция методом Лагранжа')
    plt.xlabel('X', color='gray')
    plt.ylabel('Y',color='gray')
    plt.legend(["точки", "инт. ф-ия"], loc = "best")
    return rez

def Newton(lst):
    lst = sorted(lst)
    raz = [[i[1] for i in lst]]
    for i in range(len(lst) -1):
        raz.append([])
        for j in range(1, len(lst) - i):
            raz[-1].append(raz[-2][j] - raz[-2][j-1])
    raz = [i[0] for i in raz][1:]
    h = (max([i[0] for i in lst]) - min([i[0] for i in lst])) / (len(lst) - 1)
    koef = []
    for i in range(len(lst)-1):
        if math.log10(abs(raz[i])+1) > 25:
            koef.append((int(raz[i])/fact(i+1)/h**(i+1)))
        elif raz[i] == 0:
            koef.append(0)
        else:
            koef.append((raz[i]/fact(i+1)/h**(i+1)))
    rez = []
    for x in arange(lst[0][0], lst[-1][0] + h/10, h/10):
        y = lst[0][1]
        for i in range(len(lst)-1):
            mult = 1
            for j in range(i+1):
                mult *= x - lst[j][0]
            y += koef[i] * mult
        #print(y)
        rez.append([x, y])
        

    #plt.plot([i[0] for i in lst], [i[1] for i in lst], "o")
    #plt.plot([i[0] for i in rez], [i[1] for i in rez], "r")
    
    return rez, f"{koef[-1]}(X-X0)...(X-Xn-1) + {koef[-2]}(X-X0)...(X-Xn-2) + {koef[-3]}(X-X0)...(X-Xn-3)...+ {koef[1]}(X-X0) + {koef[0]}"

def Long_Newton(lst):
    """Интерполирование методом Ньютона
            Параметры:
                    lst (list): список точек в формате [x i ,y i ]
            Возвращаемое значение: 
                    rez (list): список точек в формате [x i , y i , f i ] (fi – значение интерполированной функции в точках)
                    """
    lst2 = copy.deepcopy(lst)
    n = 15
    flag = bool(len(lst) // n)
    rez = []
    i = 0
    if flag:
        print(f"Из-за трудности в работе с столь большими данными интерполирование было выполнено с помощью кусочного интерполирования по {n} точек")
    while len(lst) > n:
        if len(lst) < n + 5:
            n = n - 5
        i += 1
        k = Newton(lst[:n])
        rez += k[0]
        lst = lst[n-1:]
    rez += Newton(lst)[0]
    plt.plot([i[0] for i in lst2], [i[1] for i in lst2], "o")
    plt.plot([i[0] for i in rez], [i[1] for i in rez], "r") 
    plt.title('Интерполяция методом Ньютона')
    plt.xlabel('X', color='gray')
    plt.ylabel('Y',color='gray')
    plt.legend(["точки", "инт. ф-ия"], loc = "best")
    return rez

def splines(lst):
    """Интерполирование методом кубических сплайнов
            Параметры:
                    lst (list): список точек в формате [x i ,y i ]
            Возвращаемое значение: 
                    rez (list): список точек в формате [x i , y i , f i ] (fi – значение интерполированной функции в точках)
                    """
    N = len(lst)
    h = [0] + [lst[i][0] - lst[i-1][0] for i in range(1, N)]
    mat = []
    mat.append([2, h[2]/(h[1]+h[2])] + [0]* (N-4) + [6*((lst[2][1] - lst[1][1])/h[2] - (lst[1][1] - lst[0][1])/h[1])/(h[2] + h[1])])
    for i in range(2, N-2):
        mat.append([0]*(N-2) + [6*((lst[i+1][1] - lst[i][1])/h[i+1] - (lst[i][1] - lst[i-1][1])/h[i])/(h[i+1] + h[i])])
        for j in range(N-1):
            if i == j:
                mat[i-1][i-1] = 2
                mat[i-1][i-2] = h[i]/(h[i]+h[i+1])
                mat[i-1][i] = h[i+1]/(h[i]+h[i+1])            
            
    mat.append([0]* (N-4) + [h[N-2]/(h[N-2]+h[N-1]), 2] + [6*((lst[-1][1] - lst[-2][1])/h[-1] - (lst[-2][1] - lst[-3][1])/h[-2])/(h[-1] + h[-2])])
    
    koef = [[0], [0], [0] + linalg.solve([i[:-1] for i in mat], [i[-1] for i in mat]).tolist() + [0], [0]]
    for i in range(1, N):
        koef[1].append(koef[2][i] * h[i] / 3 + koef[2][i-1] * h[i] / 6 + (lst[i][1] - lst[i-1][1])/h[i])
        koef[3].append((koef[2][i] - koef[2][i-1])/h[i])
        koef[0].append(lst[i][1])
    x_rez = []
    y_rez = []
    plt.plot([i[0] for i in lst], [i[1] for i in lst], "o")
    for i in range(1, N):
        step = (lst[i][0]-lst[i-1][0])/10
        x10 = []
        y10 = []
        for x in arange(lst[i-1][0], lst[i][0]+step, step):
            x10.append(x)
            y10.append(koef[3][i] * (x - lst[i][0])**3 / 6 + koef[2][i] * (x - lst[i][0])**2 / 2 + koef[1][i] * (x - lst[i][0]) + koef[0][i])
        plt.plot(x10, y10, "r")
        plt.title('Интерполяция методом кубических сплайнов')
        plt.xlabel('X', color='gray')
        plt.ylabel('Y',color='gray')
        plt.legend(["точки", "инт. ф-ия"], loc = "best")
        x_rez += x10
        y_rez += y10
    return list(zip(x_rez, y_rez))

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

def appr_gr(lst, poly, p = True, ret = False):
    summ = 0
    for i in range(len(lst)):
        lst[i].append(0)
        for j in range(len(poly)):
            lst[i][-1] += lst[i][0]**j * poly[-j-1]
        summ += (lst[i][1] - lst[i][-1])**2
    if p:
        plt.plot([i[0] for i in lst], [i[1] for i in lst], "o")
        plt.plot([i[0] for i in lst], [i[-1] for i in lst], "r")
    return [f"Дисперсия: {summ/len(lst)}", lst]

def line_appr(lst):
    """
        Линейная аппроксимация методом МНК
            Параметры:
                    lst (list): список точек в формате [x i ,y i ]
            Возвращаемое значение: 
                    None
                    """
    mat = [[0, 0, 0], [0, len(lst), 0]]
    for i in lst:
        mat[0][0] += i[0]**2
        mat[0][1] += i[0]
        mat[0][2] += i[0]*i[1]
        mat[1][0] += i[0]
        mat[1][2] += i[1]
    poly =  Jordan_dig(mat)["sol"]
    a = appr_gr(lst, poly)[0]
    plt.title('Аппроксимация линейной функцией')
    plt.xlabel('X', color='gray')
    plt.ylabel('Y',color='gray')
    plt.legend(["точки", "аппр. ф-ия"], loc = "best")
    print(f"Апроксимирующая прямая: y = {poly[0]}X+{poly[1]}")
    print(a)
    
def kvadr_appr(lst, p = True):
    """Квадратичная аппроксимация методом МНК
            Параметры:
                    lst (list): список точек в формате [x i ,y i ]
                    p (bool): технический параметр, не предполагает изменения пользователем
            Возвращаемое значение: 
                    None
                    """
    mat = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, len(lst), 0]]
    for i in lst:
        mat[0][0] += i[0]**4
        mat[0][1] += i[0]**3
        mat[0][2] += i[0]**2
        mat[0][3] += i[0]**2 * i[1]
        
        mat[1][0] += i[0]**3
        mat[1][1] += i[0]**2
        mat[1][2] += i[0]
        mat[1][3] += i[0] * i[1]
        
        mat[2][0] += i[0]**2
        mat[2][1] += i[0]
        mat[2][3] += i[1]
    
    poly = linalg.solve([i[:-1] for i in mat], [i[-1] for i in mat]).tolist()
    a = appr_gr(lst, poly, p)
    if p:
        plt.title('Аппроксимация квадратичной функцией')
        plt.xlabel('X', color='gray')
        plt.ylabel('Y',color='gray')
        plt.legend(["точки", "аппр. ф-ия"], loc = "best")
    print(f"Апроксимирующий квадратный трёхлен: y = {poly[0]}X**2+{poly[1]}X+{poly[2]}")
    print(a[0])
    if not p:
        return poly
    return a[1]

def gauss(lst):
    """Гауссова аппроксимация
            Параметры:
                    lst (list): список точек в формате [x i ,y i ]
            Возвращаемое значение: 
                    None
                    """
    B = kvadr_appr([[i[0], math.log(sys.float_info.epsilon + abs(i[1]))] for i in lst], p = False)[::-1]
    a = math.exp(B[0] - B[1]**2/4/B[2])
    b = -1/B[2]
    c = -B[1]/2/B[2]
    rez = []
    disp = 0
    for i in range(1, len(lst)):
        for j in arange(lst[i-1][0], lst[i][0], (lst[i][0]-lst[i-1][0])/10):
            y = a*math.exp(-(j - c)**2/b)
            rez.append([j, y])
        disp += (y - lst[i][1])**2
    rez.append([lst[-1][0], a*math.exp(-(lst[-1][0] - c)**2/b)])
    disp += (lst[-1][1] - a*math.exp(-(lst[-1][0] - c)**2/b))**2
    #print(rez)
    plt.plot([i[0] for i in lst], [i[1] for i in lst], "o")
    plt.plot([i[0] for i in rez], [i[-1] for i in rez], "r")
    plt.title('Аппроксимация Гауссовой функцией')
    plt.xlabel('X', color='gray')
    plt.ylabel('Y',color='gray')
    plt.legend(["точки", "аппр. ф-ия"], loc = "best")
    print(f"Апроксимирующая функция: y ={a}e^(-(x-{b})**2 / {c})")
    print("Дисперсия:", disp/len(lst))
    
    
def FFT(lst):
    """Быстрое преобразование Фурье
            Параметры:
                    lst (list): список точек в формате [x i ,y i ]
            Возвращаемое значение: 
                    None
                    """
    plt.plot([i[0] for i in lst], [i[1] for i in lst], "r")
    plt.title('Исходная функция')
    plt.xlabel('X', color='gray')
    plt.ylabel('Y',color='gray')
    plt.show()
    frequency = (len(lst)-1)/(lst[-1][0] - lst[0][0]) 
    spek = [[rfftfreq(len(lst), 1/frequency)[i], abs(rfft([j[1] for j in lst]))[i]]for i in range(1, len(lst)//2 + 1)]
    #print(spek)
    plt.plot([i[0] for i in spek], [i[1] for i in spek])    
    plt.title('Спектрограмма')
    plt.xlabel('Частота, Гц', color='gray')
    plt.ylabel('Значение',color='gray')
    plt.show()
    znach = abs(rfft([j[1] for j in lst])).tolist()
    #print(spek[1][0] - spek[0][0])
    h = spek[1][0] - spek[0][0]
    S_peak = ((znach[znach.index(max(znach))-1] + znach[znach.index(max(znach))+1])/2 + znach[znach.index(max(znach))]) * h
    S_all = - S_peak
    for i in range(1, len(znach)):
        S_all += (znach[i] + znach[i-1]) / 2 * h
    print("Коэффициент несинусоидальности", S_all/S_peak)
    
def NeDisk(lst):
    """Функция вейвлет анализа недискретными вейвлетами
            Параметры:
                    lst (list): список точек в формате [x i ,y i ]
            Возвращаемое значение: 
                    None
                    """
    lvl = 4
    fig = plt.figure(figsize=(25, 15))
    n = 2
    y = [i[1] for i in lst]
    wey = ['mexh', 'gaus1']
    for ind, wv in enumerate(wey):
        fig.add_subplot(lvl + 2, len(wey), 1 + ind).plot(y, "r")
        wavelet = pywt.ContinuousWavelet(wv, level = 1)
        coef, freqs = pywt.cwt(y, arange(1,n), wavelet)
        dDa = sum([abs(i) for i in coef[0]])
        fig.add_subplot(lvl + 2, len(wey), 3 + ind).plot(coef[0])
        plt.title(f'Вейвлет - {wv}, декомпозиция {1}-ого уровня')
        cD = 0
        #plt.show()
        #plt.matshow(coef)
        for j in range(1, lvl+1):
            plt.title(f'Вейвлет - {wv}, декомпозиция {j}-ого уровня')
            coef, freqs = pywt.cwt(coef[0], arange(1,n), wavelet)
            cD += sum([abs(j) for j in coef[0]])**2
            fig.add_subplot(lvl + 2, len(wey), (j+1)*2 + 1 + ind).plot(coef[0])
            #plt.show()
            #plt.matshow(coef)
        print(f"Коэффициент несинусоидальности вейвлетом {wv}", cD**0.5/dDa)
        
def Disk(lst):
    """Функция вейвлет анализа недискретными вейвлетами
            Параметры:
                    lst (list): список точек в формате [x i ,y i ]
            Возвращаемое значение: 
                    None
                    """
    fig = plt.figure(figsize=(25, 15))
    lvl = 4
    y = [i[1] for i in lst]
    wey = ['haar', 'db1']
    for ind, wv in enumerate(wey):
        fig.add_subplot(lvl + 2, len(wey), 1 + ind).plot(y, "r")
        plt.title('Исходная функция')
        lst2 = pywt.wavedec(y, wv, level = lvl)
        dDa = sum([abs(i) for i in lst2[0]])
        fig.add_subplot(lvl + 2, len(wey), 3 + ind).plot(lst[0])
        plt.title(f'Вейвлет - {wv}, аппроксимация')
        dD = 0
        for i in range(1, lvl+1):
            dD += sum([abs(j) for j in lst2[i]])**2
            yy = lst2[i]
            fig.add_subplot(lvl + 2, len(wey), (i+1)*2 + 1 + ind).plot(yy)
            plt.title(f'Вейвлет - {wv}, декомпозиция {i}-ого уровня')
        print(f"Коэффициент несинусоидальности вейвлетом {wv}", dD**0.5/dDa)
        
def __mama(x):
     print(x)
__mama(3)
