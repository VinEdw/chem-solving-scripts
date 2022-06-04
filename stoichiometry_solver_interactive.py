# Stoichiometry solver interactive

import chemistry_solver_base

# lists that contain the strings that translate to yes and no
yesStrs = ['yes', 'yee', 'yup', 'yep', 'y', 'yeehaw']
# noStrs = ['no', 'nay', 'nah', 'nope', 'n']

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
    print('Welcome to the stoichiometry solver.')
    print('After inputting the chemical formulas for the reactants and the products, a string representing the balanced chemical equation will be given.')
    print('Next, you will input the amount in the format "[number] [unit]" for each of the reactants/products. Leave the input blank if there is the perfect amount needed.\n')

    while True:
        print('Please name the reactants. To finish, press enter with an empty string.')
        reactants = formula_list_from_inputs()
        print('Please name the products. To finish, press enter with an empty string.')
        products = formula_list_from_inputs()
        print()

        equation = chemistry_solver_base.chemical_equation_solver(reactants, products)
        if not equation:
            print('Error: Invalid formula\n')
            continue
        
        print(chemistry_solver_base.format_chemical_equation_solution(equation))

        flipped_equation = {}
        flipped_equation['reactants'] = equation['products']
        flipped_equation['products'] = equation['reactants']

        while True:
            input_products = input('Would you like to input reactants (r) or products (p)?  ').strip().lower() in ['products', 'p']
            input_formula_list = products if input_products else reactants

            print(f'\nPlease name the amount for each of the {"products" if input_products else "reactants"} in the format "[number] [g/mol]". {"Note: Only one product amount can be inputted." if input_products else ""}')
            input_amount_list = []
            for formula in input_formula_list:
                val = input(f"{formula}  ").strip().split()
                if not bool(val):
                    continue
                if len(val) != 2:
                    print('Error: Invalid expression')
                    continue
                num, unit = val
                input_amount_list.append((float(num), unit, formula))
            if not input_amount_list:
                print('Error: You must name at least one product amount')
                continue
            print()

            if input_products:
                reactant_product_dicts = list(reversed(chemistry_solver_base.stoichiometry_solver(flipped_equation, input_amount_list)))
            else:
                reactant_product_dicts = chemistry_solver_base.stoichiometry_solver(equation, input_amount_list)
            
            if not reactant_product_dicts:
                print('Error: Invalid formula\n')
                continue

            print(chemistry_solver_base.format_stoichiometry_solver_solution(reactant_product_dicts))

            use_same_eq = input("\n\nWould you like continue using the same chemical equation?  ").strip().lower()
            if use_same_eq in yesStrs:
                continue
            else:
                break

        # ending line
        print('\n' + '-'*80 + '\n')
