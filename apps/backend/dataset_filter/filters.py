import pandas as pd
from rdkit import Chem

class Filters:

    def removeNaN(df):
        df.dropna(inplace=True)
        return df

    # -----------------------------------------------

    def removeElements(df, smilesRow):

        smiles = df[str(smilesRow)]

        valid_smiles =  'C','O','N','S','P','F','I','Br','Cl'
        flag_list = []

        for current_smiles in smiles:
            mol = Chem.MolFromSmiles(str(current_smiles))
            valid_flag = True

            for current_atom in mol.GetAtoms():
                atom = current_atom.GetSymbol()

                if atom in valid_smiles:
                    continue
                else:
                    valid_flag = False
            
            flag_list.append(valid_flag)

        df['Is_valid'] = flag_list
        df = df[(df.Is_valid == True)]
        df = df.drop(['Is_valid'], axis=1)

        return df

    # -----------------------------------------------

    def removeStrain(df, selected_organisms, organismColumn):
        
        df = df[df[str(organismColumn)].isin(list(selected_organisms))]

        return df

    # -----------------------------------------------

    def convertThreshold(df, threshold_value, threshold_unit, threshold_columns):
             
        standard_values = threshold_columns[0]
        standard_units = threshold_columns[1]
        molecular_weight = threshold_columns[2]

        df_units = list(df[standard_units])
        df_values = list(df[standard_values])
        df_MW = list(df[molecular_weight])

        new_value_list = []

        if threshold_unit == 'uM':
            for df_units_index, current_df_units in enumerate(df_units):

                if current_df_units == 'ug.mL-1':
                    new_value = (df_values[df_units_index] / df_MW[df_units_index]) * 1000

                elif current_df_units == 'uM':
                    new_value = df_values[df_units_index]

                elif current_df_units == 'nM':
                    new_value = df_values[df_units_index] * 1000
                
                new_value_list.append(new_value)

            new_units_list = ['uM'] * len(df_values)
            new_activity_list = ['Inactive'] * len(df_values)
            
            for current_value_index, current_value in enumerate(new_value_list):
                if current_value < threshold_value:
                    new_activity_list[current_value_index] = 'Active'

            df['Converted Value'] = new_value_list
            df['Converted Units'] = new_units_list
            df['Activity'] = new_activity_list

        elif threshold_unit == 'nM':
            for df_units_index, current_df_units in enumerate(df_units):

                if current_df_units == 'ug.mL-1':
                    new_value = (df_values[df_units_index] / df_MW[df_units_index])

                elif current_df_units == 'uM':
                    new_value = df_values[df_units_index] / 1000

                elif current_df_units == 'nM':
                    new_value = df_values[df_units_index]
                
                new_value_list.append(new_value)

            new_units_list = ['nM'] * len(df_values)
            new_activity_list = ['Inactive'] * len(df_values)
            
            for current_value_index, current_value in enumerate(new_value_list):
                if current_value < threshold_value:
                    new_activity_list[current_value_index] = 'Active'

            df['Converted Value'] = new_value_list
            df['Converted Units'] = new_units_list
            df['Activity'] = new_activity_list

        return df

    # -----------------------------------------------

    def checkDuplicates(df, columns):
        
        molecule_id = columns[0]
        molecule_activity = columns[1]
        molecule_value = columns[2]

        df_id = list(df[molecule_id])
        df_activity = list(df[molecule_activity])
        df_values = list(df[molecule_value])

        unique_index = []
        unique_id = []
        unique_value = []

        to_removal = []

        dup_num = df.duplicated(subset=[molecule_id]).sum()
        if dup_num != 0:

            for outer, current_id in enumerate(df_id):
                actual_index = None
                actual_id = None
                actual_value = None

                for inner in range(0, len(df_id)):
                    if outer != inner:

                        #@ Same molecule?
                        if current_id == df_id[inner]:

                            #@ Same activity?
                            if df_activity[outer] == df_activity[inner]:

                                #@ Bigger value?
                                if df_values[outer] >= df_values[inner]:
                                    actual_index = outer
                                    actual_id = current_id
                                    actual_value = df_values[outer]

                            elif df_activity[outer] != df_activity[inner]:
                                if current_id not in to_removal:
                                    to_removal.append(current_id)
                
                if (actual_index is not None):
                    try:
                        duplicate_index = unique_id.index(actual_id)

                        if actual_value > unique_value[duplicate_index]:
                            unique_value[duplicate_index] = actual_value
                            unique_index[duplicate_index] = actual_index
                        else:
                            unique_value[duplicate_index] = unique_value[duplicate_index]
                            unique_index[duplicate_index] = unique_index[duplicate_index]

                    except ValueError:
                        unique_index.append(actual_index)
                        unique_id.append(actual_id)
                        unique_value.append(actual_value)
                    
            selection = df.iloc[unique_index]
            df = df[~df[molecule_id].isin(unique_id)]
            df = df.append(selection)

            df = df[~df[molecule_id].isin(to_removal)]

        return df

    # -----------------------------------------------

    def dropColumns(df, to_drop):
        
        df = df.drop(columns=to_drop)

        return df

    # -------------------------------------------

    def shuffleRows(df):
        
        df = df.sample(frac=1)

        return df