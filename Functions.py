# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 22:34:51 2021

@author: Дима
"""

from numpy import arange
import math
import csv
import matplotlib.pyplot as plt



def srt_count(srt, dic):
    global safe_dict
    return eval(srt, dic, safe_dict)

inf_check = lambda x, k_max: x if abs(x) < k_max else "inf" 

fun = dir(math)[5:]
safe_dict = dict((k, getattr(math, k)) for k in fun)

def DIY_diff (fun_str, inter = [-1, 1], step = 0.01, filename = "derivative", graph = True):
    '''
    Выполняет численное дифференцирование на заданном интервале с заданным шагом и записываем результат в формате CSV.
    При необходимости выводит графики исходной функции и полученных значений производной
            Параметры:
                    fun_str (str): строка содержащая функцию с использованием функций доступных в библиотеке math
                    inter (list): интервал дифференцирования, список из дыух действительных чисел, второе больше или равно первого (по умолчанию [-1, 1])
                    step(float): шаг дифференцирования (по умолчанию 0.01)
                    filename(str): имя файла без расширения, куда будет сохранена таблица результатов (по умолчанию "derivative")
                    graph(bool): необходимо ли выводить график (по умолчанию True)
                    
            Возвращаемое значение:
                    [dic_lst, dic] (list(list, dic)): первый список из словарей - численное дифференцирование на интервалах без разрыва, второй словарь - все значения
    '''
    k_max = 100
    dic = {}
    dic0 = {}
    for i in arange(inter[0], inter[1] + step, step):
        check = inf_check((srt_count(fun_str, {"x":i + step}) - srt_count(fun_str, {"x":i - step})) / (2 * step), k_max)
        dic[round(i, len(str(step)) - 2)] = check
        if check == "inf":
            dic0[round(i, len(str(step)) - 2)] = "возможно, точка разрыва"
        else:
            dic0[round(i, len(str(step)) - 2)] = srt_count(fun_str, {"x":i})
        
        
    with open(filename + '.csv', 'w', newline='') as csvfile:
        fieldnames = ['X', 'Y`(X)']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=';')
        writer.writeheader()

        
        for k, v in dic.items():
            writer.writerow({
                'X':k,
                'Y`(X)':v
            })
    
    dic_lst = [{}]
    dic0_lst = [{}]
    
    for i in dic.keys():
        if dic[i] != "inf":
            dic_lst[-1][i] = dic[i]
            dic0_lst[-1][i] = dic0[i]
        elif dic[i] == "inf" and len(dic_lst[-1]) != 0:
            dic_lst.append({})
            dic0_lst.append({})

            
    if graph:      
        for i in dic0_lst:
            plt.plot(list(i.keys()), list(i.values()), "r")
        plt.title('График функции y = ' + fun_str)
        plt.xlabel('X', color='gray')
        plt.ylabel('Y',color='gray')
        plt.grid(True)
        plt.legend(['y = ' + fun_str])
        plt.show()

        for i in dic_lst:
            plt.plot(list(i.keys()), list(i.values()), "go")
        plt.title('Зависимость значения производной функции y = ' + fun_str + ' в данных точках от значения аргумента\n' )
        plt.xlabel('X', color='gray')
        plt.ylabel('Y`',color='gray')
        plt.grid(True)
        plt.legend(['y` = ' + f'({fun_str})`'])
        plt.show()
    
    
    return [dic_lst, dic]

def DIY_int(fun_str, inter = [-1, 1], step = 0.01, trapeze = True, graph = True):
    '''
    Выполняет численное интегрирование на заданном интервале с заданным шагом и записываем результат в формате CSV.
    При необходимости выводит графики исходной функции и её разбиения
            Параметры:
                    fun_str (str): строка содержащая функцию с использованием функций доступных в библиотеке math
                    inter (list): интервал интегрирования, список из дыух действительных чисел, второе больше или равно первого (по умолчанию [-1, 1])
                    step(float): шаг интегрирования (по умолчанию 0.01)
                    trapeze(bool): будет ли возможно отобразить разбиение, вне зависимости от вывода графика. Если True? детализирует вычисления в 5 раз (по умолчанию True)
                    graph(bool): необходимо ли выводить график (по умолчанию True)
                    
            Возвращаемое значение:
                    summ (float): Значение численного интеграла функции
    '''
    summ = 0
    dic = {}
    dic2 = {}
    m = 5 if trapeze else 1
    n = 0
    for i in arange(inter[0], inter[1] + step / m, step / m):
        dic[round(i, len(str(step)))] = srt_count(fun_str, {"x":i})
        if n % m == 0:
            assert DIY_diff(fun_str, [i, i], step, graph = False)[1][round(i, len(str(step)))] != "inf", \
                f"Вероятно, точка x = {round(i, len(str(step)))} является точкой разрыва, интеграл расходится" 
            dic2[round(i, len(str(step)))] = srt_count(fun_str, {"x":i})
        n += 1
        
    summ = step * (sum(list(dic2.values())) - 0.5 * float(srt_count(fun_str, {"x":inter[0]}) + srt_count(fun_str, {"x":inter[1]})))

    if graph:
        plt.plot(list(dic.keys()), list(dic.values()), "r")
        plt.plot([inter[0]]*2, [min(list(dic.values())), max(list(dic.values()))], "#30D5C8")
        plt.plot([inter[1]]*2, [min(list(dic.values())), max(list(dic.values()))], "#00A550")
        plt.plot([min(list(dic.keys())), max(list(dic.keys()))], [0]*2, "#3333FF")
        if trapeze:
            for k, v in dic2.items():
                if k != list(dic2.keys())[0]:
                    plt.plot([k]*2, [0, v], "#8B00FF")
                    plt.plot([k, k-step], [v, dic2[round(k-step, len(str(step)))]], "#8B00FF")
        plt.title('График функции y = ' + fun_str)
        plt.xlabel('X', color='gray')
        plt.ylabel('Y',color='gray')
        plt.grid(True)
        plt.legend(['y = ' + fun_str, f"Левая граница x = {inter[0]}", f"Правая граница x = {inter[1]}", "Ось OX"], fontsize=8)
        plt.show()
    return summ




def graf(lst, dic, met):
    if len(dic) == 2:
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        x = [j[-1] for j in lst]
        z = [j[0] for j in lst]
        y = [j[1] for j in lst]
        ax.plot3D(x, y, z, 'red')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
    elif len(dic) == 1:
        plt.plot([j[-1] for j in lst], [j[0] for j in lst], "o")
        plt.title(f'Решение системы {len(lst[0])-2} уравнений методом {met}')
        plt.xlabel('X', color='gray')
        plt.ylabel('Y',color='gray')
        plt.legend([i for i in dic])
    else: print("График построить невозможно")
    plt.show()

def Eiler(dic, start, inter, h = 10e-4):
    '''
    Методом Эйлера решает систему приведенных ОДУ первого порядка
            Параметры:
                    dic (dict): словарь содержащий в качестве ключей однобуквенное название функции и знак "`"  производной первого порядка\
                        а в качестве значений строка содержащая функцию с использованием функций доступных в библиотеке math
                    start (dict) словарь содержащий начальный параметр "х" - ключ и его значение. В качестве ключей однобуквенное название функции,\
                        а в качестве значений - начальные значения функции в заданной точке х
                    inter (list): интервал на котором ищутся решения, список из дыух действительных чисел, второе больше или равно первого (по умолчанию [-1, 1])
                    h (float): шаг поиска решений (по умолчанию 0.001)
                    
            Возвращаемое значение:
                    lst (list): список точек в формате [i, x i , y i , …]
    '''
    lst = [[1]]
    now = {"x":inter[0]}
    now1 = {}
    for k, v in dic.items(): 
        lst[0].append(start[k[0]] + h*srt_count(v, start))
        now[k[0]] = start[k[0]] + h*srt_count(v, start)
    for i in arange(inter[0] + h, inter[1] + h, h):
        lst.append([lst[-1][0]+1])
        for k, v in dic.items(): 
            lst[-1].append(now[k[0]] + h*srt_count(v, now))
            now1[k[0]] = now[k[0]] + h*srt_count(v, now)
        now1["x"] = round(i, 2-int(math.log10(h)))
        now = now1.copy()
        now1 = {}
        lst[-1].append(round(i, 2-int(math.log10(h))))
    graf(lst, dic, "Эйлера")
    return lst



def Eiler_koshi(dic, start, inter, h = 10e-4):
    '''
    Методом Эйлера-Коши решает систему приведенных ОДУ первого порядка
            Параметры:
                    dic (dict): словарь содержащий в качестве ключей однобуквенное название функции и знак "`"  производной первого порядка\
                        а в качестве значений строка содержащая функцию с использованием функций доступных в библиотеке math
                    start (dict) словарь содержащий начальный параметр "х" - ключ и его значение. В качестве ключей однобуквенное название функции,\
                        а в качестве значений - начальные значения функции в заданной точке х
                    inter (list): интервал на котором ищутся решения, список из дыух действительных чисел, второе больше или равно первого (по умолчанию [-1, 1])
                    h (float): шаг поиска решений (по умолчанию 0.001)
                    
            Возвращаемое значение:
                    lst (list): список точек в формате [i, x i , y i , …]
    '''
    eil = Eiler(dic, start, inter, h, dic_ret = True)
    lst = [[1]]
    now = {"x":inter[0]}
    now1 = {}
    popa = eil.pop(0)[-1]
    for k, v in dic.items(): 
        lst[0].append(start[k[0]] + h/2*(srt_count(v, start)+srt_count(v, popa)))
        now[k[0]] = start[k[0]] + h/2*(srt_count(v, start)+srt_count(v, popa))
    lst[0].append(inter[0])
    for i in arange(inter[0] + h, inter[1] + h, h):
        lst.append([lst[-1][0]+1])
        popa = eil.pop(0)[-1]
        for k, v in dic.items(): 
            lst[-1].append(now[k[0]] + h/2*(srt_count(v, now)+srt_count(v, popa)))
            now1[k[0]] = now[k[0]] + h/2*(srt_count(v, now)+srt_count(v, popa))
        
        now1["x"] = round(i, 2-int(math.log10(h)))
        now = now1.copy()
        now1 = {}
        lst[-1].append(round(i, 2-int(math.log10(h))))
    graf(lst, dic, "Эйлера-Коши")
    return lst


def Runge_kutty(dic, start, inter, h = 10e-4):
    '''
    Методом Руеге-Кутты решает систему приведенных ОДУ первого порядка
            Параметры:
                    dic (dict): словарь содержащий в качестве ключей однобуквенное название функции и знак "`"  производной первого порядка\
                        а в качестве значений строка содержащая функцию с использованием функций доступных в библиотеке math
                    start (dict) словарь содержащий начальный параметр "х" - ключ и его значение. В качестве ключей однобуквенное название функции,\
                        а в качестве значений - начальные значения функции в заданной точке х
                    inter (list): интервал на котором ищутся решения, список из дыух действительных чисел, второе больше или равно первого (по умолчанию [-1, 1])
                    h (float): шаг поиска решений (по умолчанию 0.001)
                    
            Возвращаемое значение:
                    lst (list): список точек в формате [i, x i , y i , …]
    '''
    lst = Eiler(dic, start, inter, h, dic_ret = True)
    new = lst[0][-1]
    lst = [lst[0][:-1]]
    for x in arange(inter[0] + h, inter[1] + h, h):
        lst.append([lst[-1][0]+1])
        koef = [[] for i in range(len(dic))]
        
        j = 0
        k_dic = {}
        for k, v in dic.items():
            koef[j].append(h*srt_count(v, new))
            k_dic[k] = h*srt_count(v, new)
            j += 1
            
        new1 = new.copy() 
        for i in k_dic:
            new1[i] += k_dic[i]/2
        new1["x"] += h/2
        j = 0
        for k, v in dic.items():
            koef[j].append(h*srt_count(v, new1))
            k_dic[k] = h*srt_count(v, new1)
            j += 1
        
        new1 = new.copy() 
        for i in k_dic:
            new1[i] += k_dic[i]/2
        new1["x"] += h/2
        j = 0
        for k, v in dic.items():
            koef[j].append(h*srt_count(v, new1))
            k_dic[k] = h*srt_count(v, new1)
            j += 1
            
        new1 = new.copy() 
        for i in k_dic:
            new1[i] += k_dic[i]
        new1["x"] += h
        j = 0
        for k, v in dic.items():
            koef[j].append(h*srt_count(v, new1))
            k_dic[k] = h*srt_count(v, new1)
            j += 1
            
        for l, k in enumerate(dic):
            lst[-1].append(new[k] + (1/6)*(koef[l][0] + 2*koef[l][1] + 2*koef[l][2] + koef[l][3]))
            new[k] = lst[-1][-1]
        new["x"] = round(x, 2-int(math.log10(h)))
        lst[-1].append(round(x, 2-int(math.log10(h))))
    graf(lst, dic, "Рунге-Кутты")
    return lst