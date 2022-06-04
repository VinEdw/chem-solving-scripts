# Chemistry solver base module
import math
import sympy
import periodictable as pt
import pandas as pd

def parse_formula(formula: str, index: int) -> tuple[dict[int], int]:
    formula_dict = {}
    i = index
    element = None
    element_count = 0
    # extract the elements and amounts from the chemical formula and put them in a dictionary
    while i < len(formula):
        l = formula[i]
        # if a start parentheses is encountered, call this function again on that section
        if l == '(':
            if element != None:
                formula_dict[element] = formula_dict.get(element, 0) + (element_count if element_count != 0 else 1)
                element = None
                element_count = 0
            section_dict, i = parse_formula(formula, i+1)
            for elem, amount in section_dict.items():
                formula_dict[elem] = formula_dict.get(elem, 0) + amount
            continue
        # if an end parentheses is encountered, return
        elif l == ')':
            if element != None:
                formula_dict[element] = formula_dict.get(element, 0) + (element_count if element_count != 0 else 1)
            
            if i + 1 == len(formula) or not formula[i+1].isdigit():
                i += 1
                break
            scale_factor = 0
            i += 1
            while i < len(formula) and formula[i].isdigit():
                scale_factor = scale_factor * 10 + int(formula[i])
                i += 1
            for element, amount in formula_dict.items():
                formula_dict[element] = amount * scale_factor
            break
        # if an uppercase letter is encountered, add the current element to the dictionary and start a new element
        elif l.isupper():
            if element != None:
                formula_dict[element] = formula_dict.get(element, 0) + (element_count if element_count != 0 else 1)
            element = l
            element_count = 0
        # lowercase letters continue the element started by the previous uppercase letter
        elif l.islower():
            element += l
            # if a lowercase letter does not follow another letter, return None
            if not formula[i-1].isalpha():
                return None
        # numbers tell the amount for the current element; the current element will be added to the dictionary with the amount determined so far
        elif l.isdigit():
            element_count = element_count * 10 + int(l)
        
        # if at the end of the formula, log the current element and count into the dictionary
        if i + 1 == len(formula):
            if element != None:
                formula_dict[element] = formula_dict.get(element, 0) + (element_count if element_count != 0 else 1)

        i += 1
    return (formula_dict, i)

def chemical_equation_solver(reactants: list[str], products: list[str]):
    # create dictionaries to hold the information extracted from the reactant and product lists
    # each dictionary will contain chemical formula keys that lead to a dictionary of element keys that lead to an amount for that element
    reactant_dict = {}
    product_dict = {}

    for formula_list, formula_dict in zip([reactants, products], [reactant_dict, product_dict]):
        for formula in formula_list:
            result = parse_formula(formula, 0)[0]
            if not result:
                return None
            formula_dict[formula] = result

    # these sets contain the elements in the products and reactants
    reactant_element_set = {element for elements in reactant_dict.values() for element in elements}
    product_element_set = {element for elements in product_dict.values() for element in elements}
    # if they do not match, then the chemical equation is not possible; return None
    if reactant_element_set != product_element_set:
        return None
    element_set = reactant_element_set

    # begin setting up the system of equations that will be solved to find the needed coefficients
    # Ex: two reactants & two products
    # (a_0(amount of element in formula) + a_1(amount)) / a_3 = (a_2(amount)) / a_3 + (amount for last formula)
    equation_system = []
    for element in element_set:
        equation = []
        for formula in reactants:
            num = reactant_dict[formula].get(element, 0)
            equation.append(num)
        for formula in products:
            num = product_dict[formula].get(element, 0)
            equation.append(-num)
        equation[-1] *= -1
        equation_system.append(equation)
    # solve the system of equations; it must be turned into a Matrix first to be understood by linsolve 
    equation_matrix = sympy.Matrix(equation_system)
    solution = sympy.linsolve(equation_matrix)

    # create a set of all the denominators in solution, find the lcm of the denominators, multiply each value in solution by the lcm, and stick the lcm at the end to get all the chemical formula coefficients
    denominators = {num.q for num in next(iter(solution)) if num.q != 1}
    lcm = math.lcm(*denominators) if bool(denominators) else 1
    coefficients = [int(num * lcm) for num in next(iter(solution))] + [lcm]

    # return a structure of the form {"reactants": [(amount, formula)...], "products": [(amount, formula)...]}
    return {"reactants": [(coefficients[i], formula) for i, formula in enumerate(reactants)], "products": [(coefficients[i + len(reactants)], formula) for i, formula in enumerate(products)]}

def format_chemical_equation_solution(solution) -> str:
    return f"{str.join(' + ', [((str(coefficient) + formula) if coefficient != 1 else (formula)) for coefficient, formula in solution['reactants']])} --> {str.join(' + ', [((str(coefficient) + formula) if coefficient != 1 else (formula)) for coefficient, formula in solution['products']])}"

def molar_mass_from_formula(formula: str) -> float:
    """Return the value for the molar mass of the input compound."""
    try:
        molar_mass = pt.formula(formula).mass
    except:
        return 0
    return molar_mass

def stoichiometry_solver(chemical_equation, reactants: list[tuple[float, str, str]]):
    # create a list for the reactants and products
    reactant_list = [formula for coefficient, formula in chemical_equation['reactants']]
    # product_list = [formula for coefficient, formula in chemical_equation['products']]

    # create dictionaries with helpful information on the reactants and products
    reactant_dict = {formula:{'coefficient':coefficient} for coefficient, formula in chemical_equation['reactants']}
    product_dict = {formula:{'coefficient':coefficient} for coefficient, formula in chemical_equation['products']}
    # put the molar mass of each compound into the dictionaries; if the compound is not valid, return None
    for coefficient, formula in chemical_equation['reactants']:
        molar_mass = molar_mass_from_formula(formula)
        if not molar_mass:
            return None
        reactant_dict[formula]['molar_mass'] = molar_mass
    for coefficient, formula in chemical_equation['products']:
        molar_mass = molar_mass_from_formula(formula)
        if not molar_mass:
            return None
        product_dict[formula]['molar_mass'] = molar_mass
    # if the mass/moles for a reactant was given in the function arguments, put it into the dictionary 
    for num, unit, formula in reactants:
        if formula not in reactant_list:
            return None
        if unit in ['g', 'grams', 'gram']:
            reactant_dict[formula]['mass'] = num
            reactant_dict[formula]['moles'] = num / reactant_dict[formula]['molar_mass']
        if unit in ['mol', 'moles', 'mole']:
            reactant_dict[formula]['mass'] = num * reactant_dict[formula]['molar_mass']
            reactant_dict[formula]['moles'] = num

    # find the number of reactions (times N_a) and the limiting reactant
    rxn_count, limiting_reactant = sorted([(reactant_dict[formula]['moles'] / reactant_dict[formula]['coefficient'], formula) for formula in reactant_dict if 'moles' in reactant_dict[formula]])[0]

    # use coefficient * rxn_count as the moles value for all the reactants without any moles value
    # also calculate their corresponding mass value
    for formula in reactant_dict:
        if 'moles' not in reactant_dict[formula]:
            reactant_dict[formula]['moles'] = reactant_dict[formula]['coefficient'] * rxn_count
            reactant_dict[formula]['mass'] = reactant_dict[formula]['moles'] * reactant_dict[formula]['molar_mass']

    # find the excess_mass and excess_moles for all the reactants
    for formula in reactant_dict:
        reactant_dict[formula]['excess_moles'] = reactant_dict[formula]['moles'] - reactant_dict[formula]['coefficient'] * rxn_count
        reactant_dict[formula]['excess_mass'] = reactant_dict[formula]['excess_moles'] * reactant_dict[formula]['molar_mass']

    # use coefficient*rxn_count as the moles value for all the products; also find their corresponding mass value
    for formula in product_dict:
        moles = product_dict[formula]['coefficient'] * rxn_count
        product_dict[formula]['mass'] = moles * product_dict[formula]['molar_mass']
        product_dict[formula]['moles'] = moles

    return (reactant_dict, product_dict)

def format_stoichiometry_solver_solution(solution):
    reactant_table = pd.DataFrame(solution[0]).transpose()
    reactant_table = reactant_table.astype({'coefficient':int})
    product_table = pd.DataFrame(solution[1]).transpose()
    product_table = product_table.astype({'coefficient':int})
    solution_tables_str = 'Reactants:\n' + str(reactant_table) + '\n\nProducts:\n' + str(product_table)
    return solution_tables_str

if __name__ == "__main__":
    def test_miscellaneous_working():
        a = format_chemical_equation_solution(chemical_equation_solver(['CO2', 'H2O'], ['C6H12O6', 'O2']))
        print(a)
        assert a == '6CO2 + 6H2O --> C6H12O6 + 6O2'
        b = format_chemical_equation_solution(chemical_equation_solver(['CH4', 'O2'], ['CO2', 'H2O']))
        print(b)
        assert b == 'CH4 + 2O2 --> CO2 + 2H2O'
        c = format_chemical_equation_solution(chemical_equation_solver(['N2', 'H2'], ['NH3']))
        print(c)
        assert c == 'N2 + 3H2 --> 2NH3'
        d = format_chemical_equation_solution(chemical_equation_solver(['Na', 'Cl2'], ['NaCl']))
        print(d)
        assert d == '2Na + Cl2 --> 2NaCl'
        e = format_chemical_equation_solution(chemical_equation_solver(['Cu2O', 'C'], ['Cu', 'CO2']))
        print(e)
        assert e == '2Cu2O + C --> 4Cu + CO2'
        f = format_chemical_equation_solution(chemical_equation_solver(['C2H6', 'O2'], ['CO2', 'H2O']))
        print(f)
        assert f == '2C2H6 + 7O2 --> 4CO2 + 6H2O'
        g = format_chemical_equation_solution(chemical_equation_solver(['CO', 'H2'], ['C8H18', 'H2O']))
        print(g)
        assert g == '8CO + 17H2 --> C8H18 + 8H2O'
        h = format_chemical_equation_solution(chemical_equation_solver(['Zn', 'HCl'], ['ZnCl2', 'H2']))
        print(h)
        assert h == 'Zn + 2HCl --> ZnCl2 + H2'
        i = format_chemical_equation_solution(chemical_equation_solver(['C5H8O2', 'NaH', 'HCl'], ['C5H12O2', 'NaCl']))
        print(i)
        assert i == 'C5H8O2 + 2NaH + 2HCl --> C5H12O2 + 2NaCl'
    # test_miscellaneous_working()
    def test_with_parentheses():
        a = format_chemical_equation_solution(chemical_equation_solver(['(NH4)3PO4', 'Pb(NO3)4'], ['Pb3(PO4)4', 'NH4NO3']))
        print(a)
        assert a == '4(NH4)3PO4 + 3Pb(NO3)4 --> Pb3(PO4)4 + 12NH4NO3'
        b = format_chemical_equation_solution(chemical_equation_solver(['Al2(CO3)3', 'H3PO4'], ['AlPO4', 'CO2', 'H2O']))
        print(b)
        assert b == 'Al2(CO3)3 + 2H3PO4 --> 2AlPO4 + 3CO2 + 3H2O'
        c = format_chemical_equation_solution(chemical_equation_solver(['Ca3(PO4)2', 'H2SO4'], ['CaSO4', 'Ca(H2PO4)2']))
        print(c)
        assert c == 'Ca3(PO4)2 + 2H2SO4 --> 2CaSO4 + Ca(H2PO4)2'
        d = format_chemical_equation_solution(chemical_equation_solver(['Hg(OH)2', 'H3PO4'], ['Hg3(PO4)2', 'H2O']))
        print(d)
        assert d == '3Hg(OH)2 + 2H3PO4 --> Hg3(PO4)2 + 6H2O'
    # test_with_parentheses()
    def test_three_digit_coefficients():
        a = format_chemical_equation_solution(chemical_equation_solver(['K4Fe(CN)6', 'KMnO4', 'H2SO4'], ['KHSO4', 'Fe2(SO4)3', 'MnSO4', 'HNO3', 'CO2', 'H2O']))
        print(a)
        assert a == '10K4Fe(CN)6 + 122KMnO4 + 299H2SO4 --> 162KHSO4 + 5Fe2(SO4)3 + 122MnSO4 + 60HNO3 + 60CO2 + 188H2O'
    # test_three_digit_coefficients()

    # eq_1 = {'reactants': [(16, 'HCl'), (2, 'KMnO4')], 'products': [(5, 'Cl2'), (2, 'KCl'), (2, 'MnCl2'), (8, 'H2O')]}
    # soln = stoichiometry_solver(eq_1, [(25, 'g', 'KMnO4'), (85, 'g', 'HCl')])
    # print(format_stoichiometry_solver_solution(soln))
    