import numpy as np
import uncertainties
from copy import deepcopy

def write(str):
    with open('build/ergebnis.txt', 'w') as my_file:
        my_file.write(str)
        return

def add(str):
    with open('build/ergebnis.txt', 'a') as my_file:
        my_file.write(str)
        return

def latex_float(f, dgts):
    if not type(f) == np.ndarray or type(f) == list:
        tmp = np.array([deepcopy(f)])
    else:
        tmp = deepcopy(f)
    
    res = np.array([])
    
    for elem in tmp:
        # print(type(elem))
        if type(elem) == str or type(elem) == np.str_:
            res = np.append(res, [elem])
        elif type(elem) == uncertainties.core.AffineScalarFunc or type(elem) == uncertainties.core.Variable:
            match(dgts):
                case 1:
                    float_str = "{0:.1uL}".format(elem)
                case 2:
                    float_str = "{0:.2uL}".format(elem)
                case 3:
                    float_str = "{0:.3uL}".format(elem)
                case 4:
                    float_str = "{0:.4uL}".format(elem)
                case 5:
                    float_str = "{0:.5uL}".format(elem)
                case 6:
                    float_str = "{0:.6uL}".format(elem)
                case 7:
                    float_str = "{0:.7uL}".format(elem)
                case 8:
                    float_str = "{0:.8uL}".format(elem)
                case 9:
                    float_str = "{0:.9uL}".format(elem)
                case _:
                    float_str = "{0:.2g}".format(elem)
            if r"\left(" in float_str:
                klammer, exponent = float_str.split(r"\right)")
                float_str = r"${0} \right)$ ${1}$".format(klammer, exponent)
            elif r"\pm" in float_str:
                val, unc = float_str.split(r"\pm")
                float_str = r"{0} \, \clap{{$\pm$}} \, {1}".format(val, unc)
            res = np.append(res, float_str)
        else:
            match(dgts):
                case 1:
                    float_str = "{0:.1g}".format(elem)
                case 2:
                    float_str = "{0:.2g}".format(elem)
                case 3:
                    float_str = "{0:.3g}".format(elem)
                case 4:
                    float_str = "{0:.4g}".format(elem)
                case 5:
                    float_str = "{0:.5g}".format(elem)
                case 6:
                    float_str = "{0:.6g}".format(elem)
                case 7:
                    float_str = "{0:.7g}".format(elem)
                case 8:
                    float_str = "{0:.8g}".format(elem)
                case 9:
                    float_str = "{0:.9g}".format(elem)
                case _:
                    float_str = "{0:.2g}".format(elem)
            if "e" in float_str:
                base, exponent = float_str.split("e")
                float_str =  r"{0} \clap{{$\times$}} $10^{{{1}}}$".format(base, int(exponent))
            res = np.append(res, [float_str])
    # print(res)
    return res


def tab_to_latex(arr, header, filename, cpt="Your Caption Here.", lbl="Your Label here.", dgts=2, mode="tables"):
    """
    Converts a numpy array to a LaTeX table and saves it to a file. While maintaining the originial state of the array (deepcopy).
    
    Parameters:
    arr (numpy.ndarray): The array to convert.
    header (list): The header for the table.
    filename (str): The name of the file to save the table to.
    mode: tables for normal tables or long for longtables
    dgts: int for same dgts in every column or list for column-specific
    """

    #Make a list out of arrays, because otherwise you cannot have strings in an array while also having floats.
    tmp = deepcopy(arr)
    for i in range(len(tmp)):
        if type(tmp[i]) == np.ndarray:
            tmp[i] = arr[i].tolist()
    
    arr_conv = tmp
    
    #Make Scientific Notation Compatible for LaTex-Tables
    if type(dgts) != list:
        dgts = np.ones(len(arr))*dgts
    for i in range(len(arr)):
        for j in range(len(arr[i])):
            arr_conv[i][j] = latex_float(arr[i][j], dgts[i])
    if mode == "tables":
        with open(filename, 'w') as f:
            f.write('\\begin{table}[H]\n')
            f.write('\t\\centering\n')
            f.write('\t\\caption{' + cpt + '}\n')
            f.write('\t\\label{tab:' + lbl + '}\n')
            f.write('\t\\begin{tabular}{' + ' | '.join(['c'] * len(header)) + '}\n')
            f.write('\t\t\\toprule\n')
            f.write("\t\t" + ' & '.join(header) + '\\\\\n')
            f.write('\t\t\\midrule\n\t\t')
            np.savetxt(f, np.column_stack(arr_conv), fmt = "%s", delimiter=' & ', newline=' \\\\\n\t\t ')
            f.write('\\bottomrule\n')
            f.write('\t\\end{tabular}\n')
            f.write('\\end{table}\n')
    elif mode == "long":
        with open(filename, 'w') as f:
            f.write('\\begin{longtblr}[ \n')
            f.write('\t caption = {' + cpt + '}, \n')
            f.write('\t label = {tab:' + lbl + '}, \n')
            f.write(']{ \n')
            f.write('\t colspec = {' + ' | '.join(['c'] * len(header)) + '},\n')
            f.write('\t rowhead = 1, \n')
            f.write('}\n')
            f.write('\\toprule\n')
            f.write(' & '.join(header) + '\\\\\n')
            f.write('\\midrule\n')
            np.savetxt(f, np.column_stack(arr_conv), fmt = "%s", delimiter=' & ', newline=' \\\\\n')
            f.write('\\bottomrule\n')
            f.write('\\end{longtblr}')