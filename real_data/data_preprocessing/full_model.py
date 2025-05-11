import pandas as pd
from difflib import SequenceMatcher
df = pd.read_csv('used_cars_data.csv')
columns_keep = [
    'price', 'back_legroom', 'city_fuel_economy',
    'daysonmarket',  'engine_displacement', 'front_legroom',
    'fuel_tank_volume',  'height','highway_fuel_economy', 'horsepower', 'length',
    'mileage', 'seller_rating',  'wheelbase','width', 'year',
    'body_type', 
     'listing_color',  'interior_color', 
    'fuel_type',  'engine_cylinders', 'transmission', 'wheel_system', 
     'make_name',  'owner_count', 'maximum_seating',
    'franchise_dealer', 'isCab', 'has_accidents'
]
df_filtered = df[columns_keep]
suv = df_filtered[df_filtered['body_type'] == 'SUV / Crossover']
suv.drop('body_type', axis=1, inplace=True)
suv = suv[~suv.apply(lambda x: x.astype(str).str.contains('--')).any(axis=1)]
suv = suv.dropna()
suv = suv[~suv.apply(lambda x: x.astype(str).str.contains('UNKNOWN')).any(axis=1)]
suv = suv[~suv.apply(lambda x: x.astype(str).str.contains('None')).any(axis=1)]
suv['wheel_system'] = suv['wheel_system'].replace('4WD', 'AWD')
def transform_icolor(color,unique_colors):
    color = color.upper()
    for possible_color in unique_colors:
        if possible_color in color:
            return possible_color

    return 'others'
def transform_owner_count(color):
    if color >=5:
        return 'others'
    else:
        return color
unique_colors = [color.upper() for color in suv['listing_color'].unique()]
unique_interior_color = [color.upper() for color in suv['interior_color'].unique()] 
unique_cylinders = suv['engine_cylinders'].unique()
unique_fuel_type = suv['fuel_type'].unique()
unique_transmission = suv['transmission'].unique()
unique_wheel_system = suv['wheel_system'].unique()
unique_make_name = suv['make_name'].unique()
unique_owner_count = suv['owner_count'].unique()
unique_maximum_seating = suv['maximum_seating'].unique()
suv['interior_color'] = suv['interior_color'].apply(lambda x: transform_icolor(x, unique_colors))
suv['owner_count'] = suv['owner_count'].apply(lambda x: transform_owner_count(x))
categorical_cols = suv.columns[(17-1):(28-1+1)]
suv_dummies = pd.get_dummies(suv, columns=categorical_cols)
make_name_dummies = [col for col in suv_dummies.columns if 'make_name' in col]
owner_count_dummies = [col for col in suv_dummies.columns if 'owner_count' in col]
exterior_color_dummies = [col for col in suv_dummies.columns if 'listing_color' in col]
interior_color_dummies = [col for col in suv_dummies.columns if 'interior_color' in col]
fuel_type_dummies = [col for col in suv_dummies.columns if 'fuel_type' in col]
engine_cylinders_dummies = [col for col in suv_dummies.columns if 'engine_cylinders' in col]
wheel_system_dummies = [col for col in suv_dummies.columns if 'wheel_system' in col]
transmission_dummies = [col for col in suv_dummies.columns if 'transmission' in col]
has_accidents_dummies = [col for col in suv_dummies.columns if 'has_accidents' in col]
for fuel_type_col in make_name_dummies:
    for engine_cyl_col in owner_count_dummies:
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp      
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)
for fuel_type_col in exterior_color_dummies:
    for engine_cyl_col in interior_color_dummies:
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp  
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)   
for fuel_type_col in fuel_type_dummies:
    for engine_cyl_col in engine_cylinders_dummies:
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp                                                                                    
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)              
for fuel_type_col in fuel_type_dummies:
    for wheel_system_col in wheel_system_dummies:
        interaction_col_name = f'{fuel_type_col}_X_{wheel_system_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[wheel_system_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp                                                                                     
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)    
for wheel_system_col in wheel_system_dummies:
    for engine_cyl_col in engine_cylinders_dummies:
        interaction_col_name = f'{wheel_system_col}_X_{engine_cyl_col}'
        temp = suv_dummies[wheel_system_col] * suv_dummies[engine_cyl_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp                                                                                  
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)
for wheel_system_col in transmission_dummies:
    for engine_cyl_col in engine_cylinders_dummies:
        interaction_col_name = f'{wheel_system_col}_X_{engine_cyl_col}'
        temp = suv_dummies[wheel_system_col] * suv_dummies[engine_cyl_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp  
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)           
for fuel_type_col in transmission_dummies:
    for wheel_system_col in wheel_system_dummies:
        interaction_col_name = f'{fuel_type_col}_X_{wheel_system_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[wheel_system_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp  
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)             
for fuel_type_col in fuel_type_dummies:
    for wheel_system_col in transmission_dummies:
        interaction_col_name = f'{fuel_type_col}_X_{wheel_system_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[wheel_system_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp    
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)       
for wheel_system_col in has_accidents_dummies:
    for engine_cyl_col in engine_cylinders_dummies:
        interaction_col_name = f'{wheel_system_col}_X_{engine_cyl_col}'
        temp = suv_dummies[wheel_system_col] * suv_dummies[engine_cyl_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp  
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)       
for fuel_type_col in has_accidents_dummies:
    for wheel_system_col in wheel_system_dummies:
        interaction_col_name = f'{fuel_type_col}_X_{wheel_system_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[wheel_system_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp  
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)        
for fuel_type_col in has_accidents_dummies:
    for wheel_system_col in transmission_dummies:
        interaction_col_name = f'{fuel_type_col}_X_{wheel_system_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[wheel_system_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp  
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)       
for fuel_type_col in fuel_type_dummies:
    for wheel_system_col in has_accidents_dummies:
        interaction_col_name = f'{fuel_type_col}_X_{wheel_system_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[wheel_system_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp  
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)       
for fuel_type_col in make_name_dummies:
    for engine_cyl_col in transmission_dummies:
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp  
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)
for fuel_type_col in make_name_dummies:
    for engine_cyl_col in fuel_type_dummies:
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp     
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)
for fuel_type_col in make_name_dummies:
    for engine_cyl_col in engine_cylinders_dummies:
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp       
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)
for fuel_type_col in make_name_dummies:
    for engine_cyl_col in wheel_system_dummies:
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        if all(x == 0 for x in temp) or all(x == 1 for x in temp):
            continue
        suv_dummies[interaction_col_name] =  temp     
suv_dummies.drop(suv_dummies.columns[-1], axis=1, inplace=True)
suv_dummies = suv_dummies.drop(suv_dummies.columns[[29,43,45,68,71,73,113,121,128,130,132,134]], axis=1)
df_sample = pd.DataFrame(suv_dummies)
for column in df_sample.columns[1:15]:
    df_sample[column] = df_sample[column].astype(str).str.extract(r'(\d+\.?\d*)').astype(float)
df_sample.to_csv('suvtrain.csv', index=False)