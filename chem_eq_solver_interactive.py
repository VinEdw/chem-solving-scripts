#chemical equation solver interactive

from chemistry_solver_base import chemical_equation_solver as solver
from chemistry_solver_base import format_chemical_equation_solution as format

def formula_list_from_inputs() -> list:
    i = 0
    formula_list = []
    while True:
        i += 1
        formula = input(f"{i}  ").strip()
        if not formula:
            break
        formula_list.append(formula)
    return formula_list


if __name__ == '__main__':
    print('Welcome to the chemical equation solver.')
    print('After inputting the chemical formulas for the reactants and the products, a string representing the balanced chemical equation will be given.\n')

    while True:
        print('Please name the reactants. To finish, press enter with an empty string.')
        reactants = formula_list_from_inputs()
        print('Please name the products. To finish, press enter with an empty string.')
        products = formula_list_from_inputs()
        print()

        result = solver(reactants, products)
        if not result:
            print('Error: Invalid formula\n')
            continue
        
        print(format(result) + '\n\n' + '-'*80 + '\n')
