from __future__ import division
import sys
import pandas as pd
import itertools as it
import scipy.stats as st
import collections
import lxml.etree as ET
import pickle as pk
import math
import os
import argparse
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF

########################## argparse ##########################################

def process_args(args):
    parser = argparse.ArgumentParser(usage = '%(prog)s [options]',
                                     description = 'process some value\'s'+
                                     ' genes to create a comparison\'s map.')
    parser.add_argument('-cr', '--custom_rules', 
                        type = str,
                        default = 'false',
                        choices = ['true', 'false'],
                        help = 'choose whether to use custom rules')
    parser.add_argument('-cc', '--custom_rule',
                        type = str,
                        help='custom rules to use')
    parser.add_argument('-cm', '--custom_map',
                        type = str,
                        help='custom map to use')
    parser.add_argument('-n', '--none',
                        type = str,
                        default = 'true',
                        choices = ['true', 'false'], 
                        help = 'compute Nan values')
    parser.add_argument('-pv' ,'--pValue', 
                        type = float, 
                        default = 0.1, 
                        help = 'P-Value threshold (default: %(default)s)')
    parser.add_argument('-fc', '--fChange', 
                        type = float, 
                        default = 1.5, 
                        help = 'Fold-Change threshold (default: %(default)s)')
    parser.add_argument('-td', '--tool_dir',
                        type = str,
                        required = True,
                        help = 'your tool directory')
    parser.add_argument('-op', '--option', 
                        type = str, 
                        choices = ['datasets', 'dataset_class'],
                        help='dataset or dataset and class')
    parser.add_argument('-ol', '--out_log', 
                        help = "Output log")    
    parser.add_argument('-id', '--input_data',
                        type = str,
                        help = 'input dataset')
    parser.add_argument('-ic', '--input_class', 
                        type = str, 
                        help = 'sample group specification')
    parser.add_argument('-gs', '--generate_svg',
                        type = str,
                        default = 'true',
                        choices = ['true', 'false'], 
                        help = 'generate svg map')
    parser.add_argument('-gp', '--generate_pdf',
                        type = str,
                        default = 'true',
                        choices = ['true', 'false'], 
                        help = 'generate pdf map')
    parser.add_argument('-on', '--control',
                        type = str)
    parser.add_argument('-co', '--comparison',
                        type = str, 
                        default = '1vs1',
                        choices = ['manyvsmany', 'onevsrest', 'onevsmany'])
    parser.add_argument('-ids', '--input_datas', 
                        type = str,
                        nargs = '+', 
                        help = 'input datasets')
    parser.add_argument('-na', '--names', 
                        type = str,
                        nargs = '+', 
                        help = 'input names')
    parser.add_argument('-mc',  '--choice_map',
                        type = str,
                        choices = ['HMRcoremap','ENGRO2map'])
                                    
    args = parser.parse_args()
    return args

########################### warning ###########################################

def warning(s):
    args = process_args(sys.argv)
    with open(args.out_log, 'a') as log:
            log.write(s)
            
############################ dataset input ####################################

def read_dataset(data, name):
    try:
        dataset = pd.read_csv(data, sep = '\t', header = 0, engine='python')
    except pd.errors.EmptyDataError:
        sys.exit('Execution aborted: wrong format of ' + name + '\n')
    if len(dataset.columns) < 2:
        sys.exit('Execution aborted: wrong format of ' + name + '\n')
    return dataset

############################ dataset name #####################################

def name_dataset(name_data, count):
    if str(name_data) == 'Dataset':
        return str(name_data) + '_' + str(count)
    else:
        return str(name_data)
    
############################ load id e rules ##################################

def load_id_rules(reactions):
    ids, rules = [], []
    for key, value in reactions.items():
            ids.append(key)
            rules.append(value)
    return (ids, rules)

############################ check_methods ####################################

def gene_type(l, name):
    if check_hgnc(l):
        return 'hugo_id'
    elif check_ensembl(l):
        return 'ensembl_gene_id'
    elif check_symbol(l):
        return 'symbol'
    elif check_entrez(l):
        return 'entrez_id'
    else:
        sys.exit('Execution aborted:\n' +
                 'gene ID type in ' + name + ' not supported. Supported ID'+
                 'types are: HUGO ID, Ensemble ID, HUGO symbol, Entrez ID\n')

def check_hgnc(l):
    if len(l) > 5:
        if (l.upper()).startswith('HGNC:'):
            return l[5:].isdigit()
        else:
            return False
    else:
        return False

def check_ensembl(l): 
    if len(l) == 15:
        if (l.upper()).startswith('ENS'):
            return l[4:].isdigit()
        else:  
            return False 
    else: 
        return False 

def check_symbol(l):
    if len(l) > 0:
        if l[0].isalpha() and l[1:].isalnum():
            return True
        else:
            return False
    else:
        return False

def check_entrez(l): 
    if len(l) > 0:
        return l.isdigit()
    else: 
        return False

def check_bool(b):
    if b == 'true':
        return True
    elif b == 'false':
        return False
    
############################ resolve_methods ##################################

def replace_gene_value(l, d):
    tmp = []
    err = []
    while l:
        if isinstance(l[0], list):
            tmp_rules, tmp_err = replace_gene_value(l[0], d)
            tmp.append(tmp_rules)
            err.extend(tmp_err)
        else:
            value = replace_gene(l[0], d)
            tmp.append(value)
            if value == None:
                err.append(l[0])
        l = l[1:]
    return (tmp, err)


def replace_gene(l, d):
    if l =='and' or l == 'or':
        return l
    else:
        value = d.get(l, None)
        if not(value == None or isinstance(value, (int, float))):
            sys.exit('Execution aborted: ' + value + ' value not valid\n')
        return value

def computes(val1, op, val2, cn):
    if val1 != None and val2 != None:
        if op == 'and':
            return min(val1, val2)
        else:
            return val1 + val2
    elif op == 'and':
        if cn is True:
            if val1 != None:
                return val1
            elif val2 != None:
                return val2
            else:
                return None
        else:
            return None
    else:
        if val1 != None:
            return val1
        elif val2 != None:
            return val2
        else:
            return None

def control(ris, l, cn):
    if len(l) == 1:
        if isinstance(l[0], (float, int)) or l[0] == None:
            return l[0]
        elif isinstance(l[0], list):
            return control(None, l[0], cn)
        else:
            return False
    elif len(l) > 2:
        return control_list(ris, l, cn)
    else:
        return False

def control_list(ris, l, cn):
    while l:
        if len(l) == 1:
            return False
        elif (isinstance(l[0], (float, int)) or
              l[0] == None) and l[1] in ['and', 'or']:
            if isinstance(l[2], (float, int)) or l[2] == None:
                ris = computes(l[0], l[1], l[2], cn)            
            elif isinstance(l[2], list):
                tmp = control(None, l[2], cn)
                if tmp is False:
                    return False
                else:
                    ris = computes(l[0], l[1], tmp, cn)
            else:
                return False
            l = l[3:]
        elif l[0] in ['and', 'or']:
            if isinstance(l[1], (float, int)) or l[1] == None:
                ris = computes(ris, l[0], l[1], cn)
            elif isinstance(l[1], list):
                tmp = control(None,l[1], cn)
                if tmp is False:
                    return False
                else:
                    ris = computes(ris, l[0], tmp, cn)
            else:
                return False
            l = l[2:]
        elif isinstance(l[0], list) and l[1] in ['and', 'or']:
            if isinstance(l[2], (float, int)) or l[2] == None:
                tmp = control(None, l[0], cn)
                if tmp is False:
                    return False
                else:
                    ris = computes(tmp, l[1], l[2], cn)
            elif isinstance(l[2], list):
                tmp = control(None, l[0], cn)
                tmp2 = control(None, l[2], cn)
                if tmp is False or tmp2 is False:
                    return False
                else:
                    ris = computes(tmp, l[1], tmp2, cn)
            else:
                return False
            l = l[3:]
        else:
            return False
    return ris

############################ map_methods ######################################

def fold_change(avg1, avg2):
    if avg1 == 0 and avg2 == 0:
        return 0
    elif avg1 == 0:
        return '-INF'
    elif avg2 == 0:
        return 'INF'
    else:
        return math.log(avg1 / avg2, 2)
    
def fix_style(l, col, width, dash):
    tmp = l.split(';')
    flag_col = False
    flag_width = False
    flag_dash = False
    for i in range(len(tmp)):
        if tmp[i].startswith('stroke:'):
            tmp[i] = 'stroke:' + col
            flag_col = True
        if tmp[i].startswith('stroke-width:'):
            tmp[i] = 'stroke-width:' + width
            flag_width = True
        if tmp[i].startswith('stroke-dasharray:'):
            tmp[i] = 'stroke-dasharray:' + dash
            flag_dash = True
    if not flag_col:
        tmp.append('stroke:' + col)
    if not flag_width:
        tmp.append('stroke-width:' + width)
    if not flag_dash:
        tmp.append('stroke-dasharray:' + dash)
    return ';'.join(tmp)

def fix_map(d, core_map, threshold_P_V, threshold_F_C, max_F_C):
    maxT = 12
    minT = 2
    grey = '#BEBEBE'
    blue = '#0000FF'
    red = '#E41A1C'
    for el in core_map.iter():
        el_id = str(el.get('id'))
        if el_id.startswith('R_'):
            tmp = d.get(el_id[2:])
            if tmp != None:
                p_val = tmp[0]
                f_c = tmp[1]
                if p_val < threshold_P_V:
                    if not isinstance(f_c, str):
                        if abs(f_c) < math.log(threshold_F_C, 2):
                            col = grey
                            width = str(minT)
                        else:
                            if f_c < 0:
                                col = blue
                            elif f_c > 0:
                                col = red
                            width = str(max((abs(f_c) * maxT) / max_F_C, minT))
                    else:
                        if f_c == '-INF':
                            col = blue
                        elif f_c == 'INF':
                            col = red
                        width = str(maxT)
                    dash = 'none'
                else:
                    dash = '5,5'
                    col = grey
                    width = str(minT)
                el.set('style', fix_style(el.get('style'), col, width, dash))
    return core_map

############################ make recon #######################################

def check_and_doWord(l):
    tmp = []
    tmp_genes = []
    count = 0
    while l:
        if count >= 0:
            if l[0] == '(':
                count += 1
                tmp.append(l[0])
                l.pop(0)
            elif l[0] == ')':
                count -= 1
                tmp.append(l[0])
                l.pop(0)
            elif l[0] == ' ':
                l.pop(0)
            else:
                word = []
                while l:
                    if l[0] in [' ', '(', ')']:
                        break
                    else:
                        word.append(l[0])
                        l.pop(0)
                word = ''.join(word)
                tmp.append(word)
                if not(word in ['or', 'and']):
                    tmp_genes.append(word)
        else:
            return False
    if count == 0:
        return (tmp, tmp_genes)
    else:
        return False

def brackets_to_list(l):
    tmp = []
    while l:
        if l[0] == '(':
            l.pop(0)
            tmp.append(resolve_brackets(l))
        else:
            tmp.append(l[0])
            l.pop(0)
    return tmp

def resolve_brackets(l):
    tmp = []
    while l[0] != ')':
        if l[0] == '(':
            l.pop(0)
            tmp.append(resolve_brackets(l))
        else:
            tmp.append(l[0])
            l.pop(0)
    l.pop(0)
    return tmp

def priorityAND(l):
    tmp = []
    flag = True
    while l:
        if len(l) == 1:
            if isinstance(l[0], list):
                tmp.append(priorityAND(l[0]))
            else:
                tmp.append(l[0])
            l = l[1:]
        elif l[0] == 'or':
            tmp.append(l[0])
            flag = False
            l = l[1:]
        elif l[1] == 'or':
            if isinstance(l[0], list): 
                tmp.append(priorityAND(l[0]))
            else:
                tmp.append(l[0])
            tmp.append(l[1])
            flag = False
            l = l[2:]
        elif l[1] == 'and':
            tmpAnd = []
            if isinstance(l[0], list): 
                tmpAnd.append(priorityAND(l[0]))
            else:
                tmpAnd.append(l[0])
            tmpAnd.append(l[1])
            if isinstance(l[2], list): 
                tmpAnd.append(priorityAND(l[2]))
            else:
                tmpAnd.append(l[2])
            l = l[3:]
            while l:
                if l[0] == 'and':
                    tmpAnd.append(l[0])
                    if isinstance(l[1], list): 
                        tmpAnd.append(priorityAND(l[1]))
                    else:
                        tmpAnd.append(l[1])
                    l = l[2:]
                elif l[0] == 'or':
                    flag = False
                    break
            if flag == True: #when there are only AND in list
                tmp.extend(tmpAnd)
            elif flag == False:
                tmp.append(tmpAnd)
    return tmp

def checkRule(l):
    if len(l) == 1:
        if isinstance(l[0], list):
            if checkRule(l[0]) is False:
                return False
    elif len(l) > 2:
        if checkRule2(l) is False:
            return False
    else:
        return False
    return True

def checkRule2(l):
    while l:
        if len(l) == 1:
            return False
        elif isinstance(l[0], list) and l[1] in ['and', 'or']:
            if checkRule(l[0]) is False:
                return False
            if isinstance(l[2], list):
                if checkRule(l[2]) is False:
                    return False
            l = l[3:]
        elif l[1] in ['and', 'or']:
            if isinstance(l[2], list):
                if checkRule(l[2]) is False:
                    return False
            l = l[3:]
        elif l[0] in ['and', 'or']:
            if isinstance(l[1], list):
                if checkRule(l[1]) is False:
                    return False
            l = l[2:]
        else:
            return False
    return True

def do_rules(rules):
    split_rules = []
    err_rules = []
    tmp_gene_in_rule = []
    for i in range(len(rules)):
        tmp = list(rules[i])
        if tmp:
            tmp, tmp_genes = check_and_doWord(tmp)
            tmp_gene_in_rule.extend(tmp_genes)
            if tmp is False:
                split_rules.append([])
                err_rules.append(rules[i])
            else:
                tmp = brackets_to_list(tmp)
                if checkRule(tmp):
                    split_rules.append(priorityAND(tmp))
                else:
                    split_rules.append([])
                    err_rules.append(rules[i])
        else:
            split_rules.append([])
    if err_rules:
        warning('Warning: wrong format rule in ' + str(err_rules) + '\n')
    return (split_rules, list(set(tmp_gene_in_rule)))


############################ gene #############################################

def data_gene(gene, type_gene, name, gene_custom):
    args = process_args(sys.argv)    
    for i in range(len(gene)):
        tmp = gene.iloc[i, 0]
        if tmp.startswith(' ') or tmp.endswith(' '):
            gene.iloc[i, 0] = (tmp.lstrip()).rstrip()
    gene_dup = [item for item, count in 
               collections.Counter(gene[gene.columns[0]]).items() if count > 1]
    pat_dup = [item for item, count in 
               collections.Counter(list(gene.columns)).items() if count > 1]

    if gene_dup:
        if gene_custom == None:
            if args.rules_selector == 'HMRcore':
                gene_in_rule = pk.load(open(args.tool_dir +
                                            '/local/HMRcore_genes.p', 'rb'))
            elif args.rules_selector == 'Recon':
                gene_in_rule = pk.load(open(args.tool_dir +
                                            '/local/Recon_genes.p', 'rb'))
            gene_in_rule = gene_in_rule.get(type_gene)
        else:
            gene_in_rule = gene_custom
        tmp = []
        for i in gene_dup:
            if gene_in_rule.get(i) == 'ok':
                tmp.append(i)
        if tmp:
            sys.exit('Execution aborted because gene ID '
                     +str(tmp)+' in '+name+' is duplicated\n')
    if pat_dup:
        warning('Warning: duplicated label\n' + str(pat_dup) + 'in ' + name + 
                '\n')
        
    return (gene.set_index(gene.columns[0])).to_dict()

############################ resolve ##########################################

def resolve(genes, rules, ids, resolve_none, name):
    resolve_rules = {}
    not_found = []
    flag = False
    for key, value in genes.items():
        tmp_resolve = []
        for i in range(len(rules)):
            tmp = rules[i]
            if tmp:
                tmp, err = replace_gene_value(tmp, value)
                if err:
                    not_found.extend(err)
                ris = control(None, tmp, resolve_none)
                if ris is False or ris == None:
                    tmp_resolve.append(None)
                else:
                    tmp_resolve.append(ris)
                    flag = True
            else:
                tmp_resolve.append(None)    
        resolve_rules[key] = tmp_resolve
    if flag is False:
        warning('Warning: no computable score (due to missing gene values)' +
                'for class ' + name + ', the class has been disregarded\n')
        return (None, None)
    return (resolve_rules, list(set(not_found)))

############################ split class ######################################

def split_class(classes, resolve_rules):
    class_pat = {}
    for i in range(len(classes)):
        classe = classes.iloc[i, 1]
        if not pd.isnull(classe):
            l = []
            for j in range(i, len(classes)):
                if classes.iloc[j, 1] == classe:
                    pat_id = classes.iloc[j, 0]
                    tmp = resolve_rules.get(pat_id, None)
                    if tmp != None:
                        l.append(tmp)
                    classes.iloc[j, 1] = None
            if l:
                class_pat[classe] = list(map(list, zip(*l)))
            else:
                warning('Warning: no sample found in class ' + classe +
                        ', the class has been disregarded\n')
    return class_pat

############################ map ##############################################

def maps(core_map, class_pat, ids, threshold_P_V, threshold_F_C, create_svg, create_pdf, comparison, control):
    args = process_args(sys.argv)
    if (not class_pat) or (len(class_pat.keys()) < 2):
        sys.exit('Execution aborted: classes provided for comparisons are ' +
                 'less than two\n')

    if comparison == "manyvsmany":
        for i, j in it.combinations(class_pat.keys(), 2):
            tmp = {}
            count = 0
            max_F_C = 0
            for l1, l2 in zip(class_pat.get(i), class_pat.get(j)):
                try:
                    stat_D, p_value = st.ks_2samp(l1, l2)
                    avg = fold_change(sum(l1) / len(l1), sum(l2) / len(l2))
                    if not isinstance(avg, str):
                        if max_F_C < abs(avg):
                            max_F_C = abs(avg)
                    tmp[ids[count]] = [float(p_value), avg]
                    count += 1
                except (TypeError, ZeroDivisionError):
                    count += 1
            tab = 'result/' + i + '_vs_' + j + ' (Tabular Result).tsv'
            tmp_csv = pd.DataFrame.from_dict(tmp, orient = "index")
            tmp_csv = tmp_csv.reset_index()
            header = ['ids', 'P_Value', 'Log2(fold change)']
            tmp_csv.to_csv(tab, sep = '\t', index = False, header = header)
            
            if create_svg or create_pdf:
                if args.custom_rules == 'false' or (args.custom_rules == 'true'
                                                        and args.custom_map != ''):
                    fix_map(tmp, core_map, threshold_P_V, threshold_F_C, max_F_C)
                    file_svg = 'result/' + i + '_vs_' + j + ' (SVG Map).svg'
                    with open(file_svg, 'wb') as new_map:
                        new_map.write(ET.tostring(core_map))
                        
                    
                    if create_pdf:
                        file_pdf = 'result/' + i + '_vs_' + j + ' (PDF Map).pdf'
                        renderPDF.drawToFile(svg2rlg(file_svg), file_pdf)
                    
                    if not create_svg:
                        os.remove('result/' + i + '_vs_' + j + ' (SVG Map).svg')
    elif comparison == "onevsrest":
        for single_cluster in class_pat.keys():
            t = []
            for k in class_pat.keys():
                if k != single_cluster:
                   t.append(class_pat.get(k))
            rest = []
            for i in t:
                rest = rest + i
            
            tmp = {}
            count = 0
            max_F_C = 0
            
            primo = -1

            for l1, l2 in zip(class_pat.get(single_cluster), rest):
                try:
                    stat_D, p_value = st.ks_2samp(l1, l2)
                    avg = fold_change(sum(l1) / len(l1), sum(l2) / len(l2))
                    if primo == -1:
                        primo = 0
                        print(avg)
                    if not isinstance(avg, str):
                        if max_F_C < abs(avg):
                            max_F_C = abs(avg)
                    tmp[ids[count]] = [float(p_value), avg]
                    count += 1
                except (TypeError, ZeroDivisionError):
                    count += 1
            tab = 'result/' + single_cluster + '_vs_rest (Tabular Result).tsv'
            tmp_csv = pd.DataFrame.from_dict(tmp, orient = "index")
            tmp_csv = tmp_csv.reset_index()
            header = ['ids', 'P_Value', 'Log2(fold change)']
            tmp_csv.to_csv(tab, sep = '\t', index = False, header = header)
            
            if create_svg or create_pdf:
                if args.custom_rules == 'false' or (args.custom_rules == 'true'
                                                        and args.custom_map != ''):
                    fix_map(tmp, core_map, threshold_P_V, threshold_F_C, max_F_C)
                    file_svg = 'result/' + single_cluster + '_vs_ rest (SVG Map).svg'
                    with open(file_svg, 'wb') as new_map:
                        new_map.write(ET.tostring(core_map))
                        
                    
                    if create_pdf:
                        file_pdf = 'result/' + single_cluster + '_vs_ rest (PDF Map).pdf'
                        renderPDF.drawToFile(svg2rlg(file_svg), file_pdf)
                    
                    if not create_svg:
                        os.remove('result/' + single_cluster + '_vs_ rest (SVG Map).svg') 
                        
    elif comparison == "onevsmany":
        for i, j in it.combinations(class_pat.keys(), 2):
            if i != control and j != control:
                continue
            if i == control and j == control:
                continue
            tmp = {}
            count = 0
            max_F_C = 0
            for l1, l2 in zip(class_pat.get(i), class_pat.get(j)):
                try:
                    stat_D, p_value = st.ks_2samp(l1, l2)
                    #sum(l1) da errore secondo me perchè ha null
                    avg = fold_change(sum(l1) / len(l1), sum(l2) / len(l2))
                    if not isinstance(avg, str):
                        if max_F_C < abs(avg):
                            max_F_C = abs(avg)
                    tmp[ids[count]] = [float(p_value), avg]
                    count += 1
                except (TypeError, ZeroDivisionError):
                    count += 1
            tab = 'result/' + i + '_vs_' + j + ' (Tabular Result).tsv'
            tmp_csv = pd.DataFrame.from_dict(tmp, orient = "index")
            tmp_csv = tmp_csv.reset_index()
            header = ['ids', 'P_Value', 'Log2(fold change)']
            tmp_csv.to_csv(tab, sep = '\t', index = False, header = header)
            
            if create_svg or create_pdf:
                if args.custom_rules == 'false' or (args.custom_rules == 'true'
                                                        and args.custom_map != ''):
                    fix_map(tmp, core_map, threshold_P_V, threshold_F_C, max_F_C)
                    file_svg = 'result/' + i + '_vs_' + j + ' (SVG Map).svg'
                    with open(file_svg, 'wb') as new_map:
                        new_map.write(ET.tostring(core_map))
                        
                    
                    if create_pdf:
                        file_pdf = 'result/' + i + '_vs_' + j + ' (PDF Map).pdf'
                        renderPDF.drawToFile(svg2rlg(file_svg), file_pdf)
                    
                    if not create_svg:
                        os.remove('result/' + i + '_vs_' + j + ' (SVG Map).svg')
        
        
        

    return None

############################ MAIN #############################################

def main():
    args = process_args(sys.argv)
    
    create_svg = check_bool(args.generate_svg)
    create_pdf = check_bool(args.generate_pdf)

    if os.path.isdir('result') == False:
        os.makedirs('result')

    class_pat = {}
    
    if args.option == 'datasets':
        num = 1
        for i, j in zip(args.input_datas, args.names):
            name = name_dataset(j, num)
            resolve_rules = read_dataset(i, name)
            
            resolve_rules.iloc[:, 0] = (resolve_rules.iloc[:, 0]).astype(str)
            
            ids =  pd.Series.tolist(resolve_rules.iloc[:, 0])

            resolve_rules = resolve_rules.drop(resolve_rules.columns[[0]], axis=1)
            resolve_rules = resolve_rules.replace({'None': None})
            resolve_rules = resolve_rules.to_dict('list')
            
            #Converto i valori da str a float
            to_float = lambda x: float(x) if (x != None) else None
            
            resolve_rules_float = {}
           
            for k in resolve_rules:
                resolve_rules_float[k] = list(map(to_float, resolve_rules[k])); resolve_rules_float
            
            if resolve_rules != None:
                class_pat[name] = list(map(list, zip(*resolve_rules_float.values())))
        
            num += 1
            
    if args.option == 'dataset_class':
        name = 'RAS'
        resolve_rules = read_dataset(args.input_data, name)
        resolve_rules.iloc[:, 0] = (resolve_rules.iloc[:, 0]).astype(str)
            
        ids =  pd.Series.tolist(resolve_rules.iloc[:, 0])

        resolve_rules = resolve_rules.drop(resolve_rules.columns[[0]], axis=1)
        resolve_rules = resolve_rules.replace({'None': None})
        resolve_rules = resolve_rules.to_dict('list')
            
        #Converto i valori da str a float
        to_float = lambda x: float(x) if (x != None) else None
            
        resolve_rules_float = {}
           
        for k in resolve_rules:
            resolve_rules_float[k] = list(map(to_float, resolve_rules[k])); resolve_rules_float

        classes = read_dataset(args.input_class, 'class')
        classes = classes.astype(str)
        
        if resolve_rules_float != None:
            class_pat = split_class(classes, resolve_rules_float)
    	
       
    if args.custom_rules == 'true':
        try:
            core_map = ET.parse(args.custom_map)
        except (ET.XMLSyntaxError, ET.XMLSchemaParseError):
            sys.exit('Execution aborted: custom map in wrong format')
    else:

        if args.choice_map == 'HMRcoremap':
            core_map = ET.parse(args.tool_dir+'/local/HMRcoreMap.svg')
        elif args.choice_map == 'ENGRO2map':
            core_map = ET.parse(args.tool_dir+'/local/ENGRO2map.svg')
        
    class_pat_trim = {}
    
    for key in class_pat.keys():
    	class_pat_trim[key.strip()] = class_pat[key]
        
    maps(core_map, class_pat_trim, ids, args.pValue, args.fChange, create_svg, create_pdf, args.comparison, args.control)

    print('marea alleggerito')
    print('Execution succeded')

    return None


###############################################################################

if __name__ == "__main__":
    main()
