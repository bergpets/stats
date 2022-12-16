import math
import copy
from typing import List, Union, Any

import numpy as np
import statistics
import scipy.special
import scipy.stats as st
from scipy.stats import kstest
from scipy.stats import chi2
from scipy import *
from collections import Counter
from prettytable import PrettyTable

#2 2 3 5 3 6 1 1 2 1 4 1 4 2 2 5 2 5 1 4 2 3 3 1 6 1 2 2 4 2 5 5 2 3 3 3 3 3 1 1 2 5 1 6 7 1 2 4 4 2
#3 5 8 12 20 28 30 37 39 48 50 51 55 57 58 60 61 66 67 74 76 81 88 94 106 106 107 108 108 111 111 114 116 119 125 129 141 144 146 150 152 155 157 158 161 161 167 181 187 210


def intervals_table(sample, ivs, Ni_ivs, a, b, fract_ab, h):
    ivs_glued = copy.deepcopy(ivs)
    Ni_u = copy.deepcopy(Ni_ivs)
    N = len(sample)
    n = len(ivs_glued)
    summ_Ni_u = N         
    Pi_u = []
    mi_u = []
    chi_u = []
    summ_mi_u = 0
    summ_chi_u = 0
    i = 0
    for item in ivs_glued:
        if a > sample[0] and b > sample[N - 1]:
            if  item[0] < a < item[1] and i > 0:
                item[0] = sample[0]
                for k in range(i):
                    Ni_u[i] += Ni_sample[k]
                del ivs_glued[0:i]
                del Ni_u[0:i]
            i += 1    
        elif a > sample[0] and b < sample[N - 1]:
            if item[0] < a < item[1] and i > 0:
                item[0] = sample[0]
                k = i+1
                for k in range(i):
                    Ni_u[i] += Ni_sample[k]
                del ivs_glued[0:i]
                del Ni_u[0:i]
            elif item[0] < b < item[1] and i < n:
                item[1] = sample[N - 1]
                for k in range(i, n-1):
                    Ni_u[i] += Ni_sample[k]
                del ivs_glued[i+1:len(ivs_glued)]
                del Ni_u[i+1:len(Ni_u)]
            i += 1
        
        elif a < sample[0] and b < sample[N - 1]:
            if item[0] < b < item[1] and i < n-1:
                item[1] = sample[N - 1]
                for k in range(i, n-1):
                    Ni_u[i] += Ni_sample[k]
                del ivs_glued[i+1:len(ivs_glued)]
                del Ni_u[i+1:len(Ni_u)]
            i += 1
        
        else: 
            pass

    for item in ivs_glued:
        if item[1] > a > item[0]:
            Pi_u.append(fract_ab*(item[1] - a))
        elif item[0] < b < item[1]:
            Pi_u.append(fract_ab*(b - item[0]))
        else:
            Pi_u.append(fract_ab*h)
        

    for j in range(len(ivs_glued)):
        mi_u.append(Pi_u[j]*len(sample))
        summ_mi_u += mi_u[j]

    for i in range(len(ivs_glued)):
        chi_u.append(((mi_u[i]-Ni_u[i])**2)/mi_u[i])
        summ_chi_u += chi_u[i] 
    result = [ivs_glued, Ni_u, mi_u, summ_mi_u, chi_u, summ_chi_u, Pi_u]
    return result

def get_all_div(n):
    divisors = []
    for i in range(1, int(n / 2) + 1):
        if n % i == 0:
            if i >= 7 and i <= 10:
                divisors.append(i)
    return divisors


def inter_input() -> list[Union[list[list[float]], list[int]]]:
    lst = []
    ni = []
    n = int(input("введите количество интервалов: "))
    for i in range(n):
        interval = input("Интервал: ")
        interval = sorted(list(map(float, interval.split())))
        lst.append(interval)
    print()
    for i in range(n):
        Ni = int(input("Ni: "))
        ni.append(Ni)

    res = [lst, ni]
    return res


def get_all_div(n) -> list[int]:
    divisors = []
    for i in range(1, int(n / 2) + 1):
        if n % i == 0:
            if 7 <= i <= 10:
                divisors.append(i)
    return divisors


if __name__ == "__main__":
    # 3 5 8 12 20 28 30 37 39 48 50 51 55 57 58 60 61 66 67 74 76 81 88 94 106 106 107 108 108 111 111 114 116 119 125 129 141 144 146 150 152 155 157 158 161 161 167 181 187 210
    # 0 8 4 2 3 3 10 1 9 1 1 4 14 0 11 6 1 2 11 2 0 3 5 3 12 2 11 1 4 0 4 0 14 2 10 2 5 0 4 1 0 2 4 10 6 3 1 5 3 2
    # 370 400 410 410 360 430 440 430 430 430 430 420 420 410 400 390 400 400 380 440 440 450 440 380 420 420 400 400 410 370 390 390 430 420 410 390 440 420 430 410 400 390 370 380 360 380 360 340 320 300

    filename = input("Имя файла: ")
    f = open(filename, "w", encoding="UTF-8")

    # вводим выборку и сортируем ее по возрастанию
    sample = input("Введите вашу выборку строкой (через пробел): ")
    sample = sorted(list(map(float, sample.split())))

    # выборка уникальных элементов, частота и Ni
    un_sample = sorted(list(set(sample)))
    sample_and_n = []
    sample_and = Counter(sample)
    sample_and = sorted(sample_and.items())
    for item in sample_and:
        sample_and_n.append(list(item))

    Ni_sample = []
    for item in sample_and_n:
        i = item[1]
        Ni_sample.append(i)

    # относительная частота вариации
    wi = []
    for i in range(len(un_sample)):
        wi.append(float(Ni_sample[i]) / float(len(sample)))

    tb1 = PrettyTable()
    tb1.field_names = ["Xi", "Ni", "Wi"]
    for i in range(len(un_sample)):
        tb1.add_row([un_sample[i], Ni_sample[i], wi[i]])

    f.write(str(tb1))
    print(tb1)

    # среднее выборочное, дисперсия, сигма выборки, s^2, s
    x = 0
    for i in sample:
        x += i
    x = x / len(sample)
    f.write(f"\nCреднее выборочное Хв: {x}\n")

    D = 0
    for i in range(len(un_sample)):
        D += (float(un_sample[i]) ** 2) * float(Ni_sample[i])
    D = D / len(sample) - x**2

    sigm = D ** (1 / 2)
    f.write(f"Dв: {round(D, 2)}\n")
    f.write(f"Выборочное среднеквадратическое отклонение sigm_в: {round(sigm, 2)}\n")

    s_2 = round(len(sample) / (len(sample) - 1) * D, 2)
    s = s_2 ** (1 / 2)
    f.write(f"S^2: {s_2}\n")
    f.write(f"S: {round((s_2) ** (1 / 2), 2)}\n")

    # функция распределения
    f.write("\nФункция распределения F(x):\n\n")
    f.write(f"0\tx <= {un_sample[0]}\n")
    f_distr = [0]
    for i in range(len(un_sample) - 1):
        f_distr[i] += wi[i]
        a = f_distr[i]
        f_distr.append(a)
        f.write(f"{round(f_distr[i], 2)}\t{un_sample[i]} < x <= {un_sample[i + 1]}\n")
    f.write(f"1\t{un_sample[len(un_sample) - 1]} < x\n")

    # доверительные интервалы при известном и неизвестном мат ожидании
    # доверительные интервалы при известном и неизвестном а
    # P=[round(v,1) for v in P]
    f.write("\nДОВЕРИТЕЛЬНЫЕ ИНТЕРВАЛЫ для 0.9 и 0.95\n")
    f.write("\nДоверительные интервалы при известном мат ожидании:\n")
    interval_09 = list(
        st.norm.interval(confidence=0.9, loc=np.mean(sample), scale=st.sem(sample))
    )
    interval_09 = [round(v, 2) for v in interval_09]
    f.write(f"При gam = 0.9: {interval_09}\n")
    interval_095 = list(
        st.norm.interval(confidence=0.95, loc=np.mean(sample), scale=st.sem(sample))
    )
    interval_095 = [round(v, 2) for v in interval_095]
    f.write(f"При gam = 0.95: {interval_095}\n")

    f.write("\nДоверительные интервалы при неизвестном мат ожидании:\n")
    interval_09 = list(st.norm.interval(confidence=0.9, loc=np.mean(sample)))
    interval_09 = [round(v, 2) for v in interval_09]
    f.write(f"При gam = 0.9: {interval_09}\n")
    interval_095 = list(st.norm.interval(confidence=0.95, loc=np.mean(sample)))
    interval_095 = [round(v, 2) for v in interval_095]
    f.write(f"При gam = 0.95: {interval_095}\n")

    f.write("\nДоверительные интервалы считая а известным:\n")
    f.write(f"a = [xв] = {round(x, 0)}\n")

    x_left_09 = chi2.ppf(1 - (1 - 0.9) / 2, len(sample))
    x_right_09 = chi2.ppf(1 - (1 + 0.9) / 2, len(sample))
    x_left_095 = chi2.ppf(1 - (1 - 0.95) / 2, len(sample))
    x_right_095 = chi2.ppf(1 - (1 + 0.95) / 2, len(sample))

    summ_a = 0
    for i in sample:
        summ_a += (i - round(x, 0)) ** 2

    f.write(
        f"При gam = 0.9: [{round((summ_a / x_left_09) ** (1 / 2), 2)}, {round((summ_a / x_right_09) ** (1 / 2), 2)}]\n"
    )
    f.write(
        f"При gam = 0.95: [{round((summ_a / x_left_095) ** (1 / 2), 2)}, {round((summ_a / x_right_095) ** (1 / 2), 2)}]\n"
    )

    f.write("\nДоверительные интервалы считая а неизвестным:\n\n")

    x1_left_09 = chi2.ppf(1 - (1 - 0.9) / 2, len(sample) - 1)
    x1_right_09 = chi2.ppf(1 - (1 + 0.9) / 2, len(sample) - 1)
    x1_left_095 = chi2.ppf(1 - (1 - 0.95) / 2, len(sample) - 1)
    x1_right_095 = chi2.ppf(1 - (1 + 0.95) / 2, len(sample) - 1)

    left_09 = ((len(sample) - 1) * s_2) / x1_left_09
    right_09 = ((len(sample) - 1) * s_2) / x1_right_09
    left_095 = ((len(sample) - 1) * s_2) / x1_left_095
    right_095 = ((len(sample) - 1) * s_2) / x1_right_095

    f.write(
        f"При gam = 0.9: [{round((left_09) ** (1 / 2), 2)}, {round((right_09) ** (1 / 2), 2)}]\n"
    )
    f.write(
        f"При gam = 0.95: [{round((left_095) ** (1 / 2), 2)}, {round((right_095) ** (1 / 2), 2)}]\n"
    )

    # интервальный вариационный ряд
    f.write("\nИнтервальный вариационный ряд\n")
    ivs_n = get_all_div(int(sample[len(sample) - 1] - sample[0]))
    print(f"\nРекомендуемое количество интервалов: {ivs_n}")
    n = int(input("Введите желаемое количество интервалов для ИВР: "))
    h = (sample[len(sample) - 1] - sample[0]) / n
    ivs = []
    i = 0
    while i != n:
        segment = [un_sample[0] + h * i, un_sample[0] + h * (i + 1)]
        ivs.append(segment)
        i += 1
    f.write(f"\nШаг интервала h: {h}\n")

    avg_ivs = []
    for item in ivs:
        i = 0
        avg_ivs.append(item[i] + (item[1] - item[0]) / 2)
        i += 1

    Ni_ivs = []
    for i in range(n + 1):
        Ni_ivs.append(0)
    i = 0
    for coup in ivs:
        for x_ in sample_and_n:
            if x_[0] > coup[0] and x_[0] < coup[1]:
                Ni_ivs[i] += x_[1]
            elif x_[0] == coup[0] and i == 0:
                Ni_ivs[i] += x_[1]
            elif x_[0] == coup[0] and i != 0:
                pass
            elif x_[0] == coup[1] and i < len(ivs) and x_[1] % 2 == 0:
                Ni_ivs[i + 1] += x_[1] / 2
                Ni_ivs[i] += x_[1] / 2
            elif x_[0] == coup[1] and i < len(ivs) and x_[1] % 2 == 1:
                Ni_ivs[i + 1] += x_[1] / 2 + 0.5
                Ni_ivs[i] += x_[1] / 2 - 0.5
            elif x_[0] == coup[1] and i == len(ivs):
                Ni_ivs[i] += x_[1]
        i += 1
    Ni_ivs[len(Ni_ivs) - 2] += Ni_ivs[len(Ni_ivs) - 1]
    Ni_ivs.pop()

    wi_norm = []
    for item in Ni_ivs:
        wi_norm.append(round(item / len(sample), 3))

    tb2 = PrettyTable()
    tb2.field_names = ["Ii", "Xc", "Ni", "Wi"]
    for i in range(len(avg_ivs)):
        tb2.add_row([ivs[i], avg_ivs[i], Ni_ivs[i], wi_norm[i]])

    f.write(str(tb2))

    x_ivs = 0
    i = 0
    for item in avg_ivs:
        x_ivs += item * Ni_ivs[i]
        i += 1
    x_ivs = x_ivs / len(sample)

    D_ivs = 0
    for i in range(len(avg_ivs)):
        D_ivs += (float(avg_ivs[i]) ** 2) * float(Ni_ivs[i])
    D_ivs = D_ivs / len(sample) - x_ivs**2
    sigm_ivs = D_ivs ** (1 / 2)

    f.write(f"\nСреднее выборочное ИВР x: {x_ivs}\n")
    f.write(f"Дисперсия ИВР D: {D_ivs}\n")
    f.write(f"Сигма sigm: {sigm_ivs}\n")

    # ассиметрия и эксцесс
    asm = 0
    i = 0
    for item in avg_ivs:
        asm += ((item - x_ivs) ** 3) * Ni_ivs[i]
        i += 1
    asm = asm / (len(sample) * sigm_ivs**3)
    exc = 0
    i = 0
    for item in avg_ivs:
        exc += ((item - x_ivs) ** 4) * Ni_ivs[i]
        i += 1
    exc = exc / (len(sample) * sigm_ivs**4) - 3
    f.write(f"\nАссиметрия и эксцесс\n")
    f.write(f"Ассиметрия a: {round(asm, 2)}\n")
    f.write(f"Эксцесс e: {round(exc, 2)}\n")

    # выравнивание по нормальному закону
    mi_norm = []
    Pi_norm = []
    chi_norm = []
    summ_mi_n = 0
    summ_Ni_n = 0
    summ_chi_n = 0
    i = 0

    for item in ivs:
        Pi_norm.append(
            st.norm.cdf((item[1] - x_ivs) / sigm_ivs)
            - st.norm.cdf((item[0] - x_ivs) / sigm_ivs)
        )
        i += 1

    for j in range(len(ivs)):
        mi_norm.append(Pi_norm[j] * len(sample))
        summ_mi_n += mi_norm[j]
        summ_Ni_n += Ni_ivs[j]

    for i in range(len(ivs)):
        chi_norm.append(((mi_norm[i] - Ni_ivs[i]) ** 2) / mi_norm[i])
        summ_chi_n += chi_norm[i]
    f.write("\nВыравнивание частот по нормальному закону:\n")
    tb3 = PrettyTable()
    tb3.field_names = ["Ii", "Ni", "mi", "chi^2"]
    for i in range(len(ivs)):
        tb3.add_row([ivs[i], Ni_ivs[i], mi_norm[i], chi_norm[i]])
    tb3.add_row(["summ:", summ_Ni_n, round(summ_mi_n, 2), summ_chi_n])

    f.write(str(tb3))
    f.write(f"\nPi: {Pi_norm}\n")

    f.write("\nВыравнивание частот по равномерному закону (По Пирсону): \n")
    r1 = int(2 * x * 100)
    r2 = int(round(3.46 * sigm, 2) * 100)

    L = np.array([[100, 100], [-100, 100]])
    R = np.array([r1, r2])
    X1 = np.linalg.inv(L).dot(R)
    X1 = list([round(i, 2) for i in X1])
    int_a = X1[0]
    int_b = X1[1]
    f.write(f"Отрезок (a;b) равен: {X1}\n")
    fract_ab = 1 / (int_b - int_a)
    f.write(f"Коэффициент 1/(b-a): {fract_ab}\n")

    res = intervals_table(sample, ivs, Ni_ivs, int_a, int_b, fract_ab, h)

    ivs_glued = res[0]
    Ni_unif = res[1]
    mi_unif = res[2]
    summ_mi_u = res[3]
    chi_unif = res[4]
    summ_chi_u = res[5]
    Pi_unif = res[6]

    wi_unif = []
    for item in Ni_unif:
        wi_unif.append(item / len(sample))

    tb4 = PrettyTable()
    tb4.field_names = ["Ii", "Ni", "mi", "chi^2"]
    for i in range(len(ivs_glued)):
        tb4.add_row([ivs_glued[i], Ni_unif[i], mi_unif[i], chi_unif[i]])
    tb4.add_row(["summ:", summ_Ni_n, round(summ_mi_u, 2), summ_chi_u])

    f.write(str(tb4))
    f.write(f"\nPi: {Pi_unif}\n")

    # выравниваение частот по методу максимального правдоподобия

    f.write(
        "\nВыравниваение частот по методу максимального правдоподобия: \n",
    )
    mmp_a = sample[0]
    mmp_b = sample[len(sample) - 1]
    f.write(f"Отрезок (a;b) равен: [{mmp_a}, {mmp_b}]\n")
    fract_ab_mmp = 1 / (mmp_b - mmp_a)
    f.write(f"Коэффициент 1/(b-a): {fract_ab_mmp}\n")

    Pi_mmp = []
    mi_mmp = []
    Ni_mmp = copy.deepcopy(Ni_ivs)
    chi_mmp = []

    summ_mi_m = 0
    summ_Ni_m = len(sample)
    summ_chi_m = 0

    for i in range(len(ivs)):
        Pi_mmp.append(fract_ab_mmp * h)
        mi_mmp.append(Pi_mmp[i] * len(sample))
        summ_mi_m += mi_mmp[i]
        chi_mmp.append(((mi_mmp[i] - Ni_mmp[i]) ** 2) / mi_mmp[i])
        summ_chi_m += chi_mmp[i]
        i += 1

    tb5 = PrettyTable()
    tb5.field_names = ["Ii", "Ni", "mi", "chi^2"]
    for i in range(len(ivs)):
        tb5.add_row([ivs[i], Ni_mmp[i], mi_mmp[i], chi_mmp[i]])
    tb5.add_row(["summ:", summ_Ni_m, round(summ_mi_m, 2), summ_chi_m])
    f.write(str(tb5))
    f.write(f"\nPi: {Pi_mmp}\n")

    # несмещенная оценка с минимальной дисперсией

    f.write("\nНесмещенная оценка с минимальной дисперсией\n")
    min_a = round(
        (len(sample) * sample[0] - sample[len(sample) - 1]) / (len(sample) - 1), 2
    )
    min_b = round(
        (len(sample) * sample[len(sample) - 1] - sample[0]) / (len(sample) - 1), 2
    )
    f.write(f"Отрезок (a;b) равен: [{min_a}, {min_b}]\n")
    fract_ab_min = 1 / (min_b - min_a)
    f.write(f"Коэффициент 1/(b-a): {fract_ab_min}\n")

    res1 = intervals_table(sample, ivs, Ni_ivs, min_a, min_b, fract_ab_min, h)

    ivs_glued_m = res1[0]
    print(ivs_glued_m)
    Ni_min = res1[1]
    mi_min = res1[2]
    summ_mi_mn = res1[3]
    chi_min = res1[4]
    summ_chi_mn = res1[5]
    Pi_min = res1[6]

    tb6 = PrettyTable()
    tb6.field_names = ["Ii", "Ni", "mi", "chi^2"]
    for i in range(len(ivs_glued_m)):
        tb6.add_row([ivs_glued_m[i], Ni_min[i], mi_min[i], chi_min[i]])
    tb6.add_row(["summ:", len(sample), round(summ_mi_mn, 2), summ_chi_mn])
    f.write(str(tb6))
    f.write(f"\nPi: {Pi_min}\n")

    # таблица значений по равномерному закону:

    f.write("\n Таблица значений по равномерному закону: \n")

    tb7 = PrettyTable()
    tb7.field_names = ["Параметр", "М. Пирсона", "ММП", "Несмещ-ая оценка с мин. дисп."]
    tb7.add_row(["a", int_a, mmp_a, min_a])
    tb7.add_row(["b", int_b, mmp_b, min_b])
    tb7.add_row(["1/(b-a)", fract_ab, fract_ab_mmp, fract_ab_min])
    tb7.add_row(["chi^2", summ_chi_u, summ_chi_m, summ_chi_mn])
    f.write(str(tb7))

    # Критерий Колмогорова

    f.write("\nКритерий Колмогорова для нормального распределения: \n")

    f_ivs = [0]
    f_norm = [0]
    diff = []
    for i in range(len(ivs)):
        f_norm[i] += round(Pi_norm[i], 2)
        p = f_norm[i]
        f_norm.append(p)
    f_norm.pop()

    for i in range(len(ivs)):
        f_ivs[i] += wi_norm[i]
        a = f_ivs[i]
        f_ivs.append(a)
    f_ivs.pop()

    for i in range(len(ivs)):
        diff.append(abs(f_ivs[i] - f_norm[i]))

    tb8 = PrettyTable()
    tb8.field_names = ["F(x)", "Pi_norm", "Fnorm", "|F(x) - Fnorm|"]
    for i in range(len(ivs)):
        tb8.add_row(
            [
                round(f_ivs[i], 2),
                round(Pi_norm[i], 2),
                round(f_norm[i], 2),
                round(diff[i], 2),
            ]
        )
    f.write(str(tb8))

    f.write(f"\nDnorm = max(|F(x) - Fnorm|): {max(diff)}\n")
    f.write(f"lambda_norm = Dnorm * N^0.5 : {max(diff) * ((len(sample)) ** (0.5))}\n")
    kstnorm = scipy.special.kolmogorov(max(diff) * ((len(sample)) ** (0.5)))
    if 0.9 < kstnorm < 1:
        f.write(f"Распределение может быть рассмотрено: P = {kstnorm}\n")
    else:
        f.write(
            f"Распределение не может быть рассмотрено: P = {kstnorm}\n",
        )

    f.write(
        "\nКритерий Колмогорова для равн-ого распределения по Методу Пирсона: \n\n",
    )

    f_unif = [0]
    f_n = [0]
    diff1 = []
    for i in range(len(ivs_glued)):
        f_unif[i] += round(Pi_unif[i], 2)
        u = f_unif[i]
        f_unif.append(u)
    f_unif.pop()

    for i in range(len(ivs_glued)):
        f_n[i] += wi_unif[i]
        w = f_n[i]
        f_n.append(w)
    f_n.pop()

    for i in range(len(ivs_glued)):
        diff1.append(abs(f_n[i] - f_unif[i]))

    tb9 = PrettyTable()
    tb9.field_names = ["F(x)", "Pi_ravn", "Fravn", "|F(x) - Fravn|"]
    for i in range(len(ivs_glued)):
        tb9.add_row(
            [
                round(f_n[i], 2),
                round(Pi_unif[i], 2),
                round(f_unif[i], 2),
                round(diff1[i], 2),
            ]
        )
    f.write(str(tb9))

    f.write(f"\nDravn = max(|F(x) - Fravn|): {max(diff1)}\n")
    f.write(
        f"lambda_ravn = Dravn * N^0.5 : {max(diff1) * ((len(sample)) ** (0.5))}\n",
    )

    kstunif = scipy.special.kolmogorov(max(diff1) * ((len(sample)) ** (0.5)))
    if 0.9 < kstunif < 1:
        f.write(f"Распределение может быть рассмотрено: P = {kstunif}\n")
    else:
        f.write(
            f"Распределение не может быть рассмотрено: P = {kstunif}\n",
        )

    # подсчет степеней свободы и уровней значимости

    f.write("\nПодсчет степеней свободы и уровней значимости: \n")

    df_norm = len(ivs) - 2 - 1
    df_unif = len(ivs_glued) - 2 - 1
    df_mmp = len(ivs) - 2 - 1
    df_min = len(ivs_glued_m) - 2 - 1

    f.write(f"Степень свободы для нормального распределния: {df_norm}\n")
    f.write(f"Степень свободы для равномерного по Пирсону: {df_unif}\n")
    f.write(f"Степень свободы для равномерного по ММП: {df_mmp}\n")
    f.write(
        f"Степень свободы для равномерного по несмещ. оценке с мин. дисп.: {df_min}\n"
    )
    f.write("\n")

    if chi2.ppf(0.95, df_norm) >= summ_chi_n:
        f.write(
            f"Гипотезу о нормальном распределении можно принять: {chi2.ppf(0.95, df_norm)} > {summ_chi_n}\n"
        )
    else:
        f.write(
            f"Гипотеза о нормальном не принимается: {chi2.ppf(0.95, df_norm)} < {summ_chi_n}\n"
        )

    if chi2.ppf(0.95, df_unif) >= summ_chi_u:
        f.write(
            f"Гипотезу о равномерном по Пирсону можно принять: {chi2.ppf(0.95, df_unif)} > {summ_chi_u}\n"
        )
    else:
        f.write(
            f"Гипотеза o равномерном по Пирсону не принимается: {chi2.ppf(0.95, df_unif)} < {summ_chi_u}\n"
        )

    if chi2.ppf(0.95, df_mmp) >= summ_chi_m:
        f.write(
            f"Гипотезу о равномерном по ММП можно принять: {chi2.ppf(0.95, df_mmp)} > {summ_chi_m}\n"
        )
    else:
        f.write(
            f"Гипотеза о равномерном по ММП не принимается: {chi2.ppf(0.95, df_mmp)} < {summ_chi_m}\n"
        )

    if chi2.ppf(0.95, df_min) >= summ_chi_mn:
        f.write(
            f"Гипотезу о равномерном по несмещ. оценке с мин. дисп. можно принять: {chi2.ppf(0.95, df_min)} > {summ_chi_mn}\n"
        )
    else:
        f.write(
            f"Гипотеза о равномерном по несмещ. оценке с мин. дисп. не принимается: {chi2.ppf(0.95, df_min)} < {summ_chi_mn}\n"
        )

    f.close()
    print(f"Файл '{filename}' готов")
