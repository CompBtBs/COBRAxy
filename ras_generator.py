from __future__ import division
import sys
import pandas as pd
import collections
import pickle as pk
import math
import argparse

########################## argparse ##########################################

def process_args(args):
    parser = argparse.ArgumentParser(usage = '%(prog)s [options]',
                                     description = 'process some value\'s'+
                                     ' genes to create a comparison\'s map.')
    parser.add_argument('-rs', '--rules_selector', 
                        type = str,
                        default = 'HMRcore',
                        choices = ['HMRcore', 'Recon', 'ENGRO2','Custom'], 
                        help = 'chose which type of dataset you want use')
    parser.add_argument('-cr', '--custom',
                        type = str,
                        help='your dataset if you want custom rules')
    parser.add_argument('-n', '--none',
                        type = str,
                        default = 'true',
                        choices = ['true', 'false'], 
                        help = 'compute Nan values')
    parser.add_argument('-td', '--tool_dir',
                        type = str,
                        required = True,
                        help = 'your tool directory')
    parser.add_argument('-ol', '--out_log', 
                        help = "Output log")    
    parser.add_argument('-id', '--input',
                        type = str,
                        help = 'input dataset')
    parser.add_argument('-ra', '--ras_output',
                        type = str,
                        required = True,
                        help = 'ras output')
    
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

def make_recon(data):
    try:
        import cobra as cb
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            recon = cb.io.read_sbml_model(data)
        react = recon.reactions
        rules = [react[i].gene_reaction_rule for i in range(len(react))]
        ids = [react[i].id for i in range(len(react))]
    except:
        try:
            data = (pd.read_csv(data, sep = '\t', dtype = str, engine='python')).fillna('')
            if len(data.columns) < 2:
                sys.exit('Execution aborted: wrong format of '+
                         'custom datarules\n')
            if not len(data.columns) == 2:
                warning('Warning: more than 2 columns in custom datarules.\n' +
                        'Extra columns have been disregarded\n')
            ids = list(data.iloc[:, 0])
            rules = list(data.iloc[:, 1])
        except pd.errors.EmptyDataError:
            sys.exit('Execution aborted: wrong format of custom datarules\n')
        except pd.errors.ParserError:
            sys.exit('Execution aborted: wrong format of custom datarules\n')            
    split_rules, tmp_genes = do_rules(rules)
    gene_in_rule = {}
    for i in tmp_genes:
        gene_in_rule[i] = 'ok'
    return (ids, split_rules, gene_in_rule)

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
            elif args.rules_selector == 'ENGRO2':
                gene_in_rule = pk.load(open(args.tool_dir +
                                            '/local/ENGRO2_genes.p', 'rb'))
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

############################ create_ras #######################################

def create_ras (resolve_rules, dataset_name, rules, ids, file):

    if resolve_rules == None:
        warning("Couldn't generate RAS for current dataset: " + dataset_name)

    for geni in resolve_rules.values():
        for i, valori in enumerate(geni):
            if valori == None:
                geni[i] = 'None'
                
    output_ras = pd.DataFrame.from_dict(resolve_rules)
    
    output_ras.insert(0, 'Reactions', ids)
    output_to_csv = pd.DataFrame.to_csv(output_ras, sep = '\t', index = False)
    
    text_file = open(file, "w")
    
    text_file.write(output_to_csv)
    text_file.close()

############################ MAIN #############################################

def main():
    args = process_args(sys.argv)

    if args.rules_selector == 'HMRcore':        
        recon = pk.load(open(args.tool_dir + '/local/HMRcore_rules.p', 'rb'))
    elif args.rules_selector == 'Recon':
        recon = pk.load(open(args.tool_dir + '/local/Recon_rules.p', 'rb'))
    elif args.rules_selector == 'ENGRO2':
        recon = pk.load(open(args.tool_dir + '/local/ENGRO2_rules.p', 'rb'))
    elif args.rules_selector == 'Custom':
        ids, rules, gene_in_rule = make_recon(args.custom)
        
    resolve_none = check_bool(args.none)
    
    
    name = "RAS Dataset"
    dataset = read_dataset(args.input, "dataset")

    dataset.iloc[:, 0] = (dataset.iloc[:, 0]).astype(str)

    type_gene = gene_type(dataset.iloc[0, 0], name) 
        
    if args.rules_selector != 'Custom':
        genes = data_gene(dataset, type_gene, name, None)
        ids, rules = load_id_rules(recon.get(type_gene))
    elif args.rules_selector == 'Custom':
        genes = data_gene(dataset, type_gene, name, gene_in_rule)
    
    resolve_rules, err = resolve(genes, rules, ids, resolve_none, name)

    create_ras(resolve_rules, name, rules, ids, args.ras_output)
      
    if err != None and err:
        warning('Warning: gene\n' + str(err) + '\nnot found in class '
            + name + ', the expression level for this gene ' +
            'will be considered NaN\n')

    
    print('Execution succeded')

    return None

###############################################################################

if __name__ == "__main__":
    main()
