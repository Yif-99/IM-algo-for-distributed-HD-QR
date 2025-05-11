import pandas as pd
from difflib import SequenceMatcher
import json
df = pd.read_csv('used_cars_data.csv')
columns_keep = [
    'price', 'back_legroom', 'city_fuel_economy',
    'daysonmarket',  'engine_displacement', 'front_legroom',
    'fuel_tank_volume',  'height','highway_fuel_economy', 'horsepower', 'length',
    'mileage', 'seller_rating',  'wheelbase','width', 'year',
    'body_type', 'listing_color',  'interior_color', 
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
def transform_color(color):
    color = color.upper() 
    return color
def transform_engine_cylinders(aa):
    aa = aa.upper()  
    unique_cylinders = ['I3','I4','I5','I6','V6','V8','H4','H6']
    for possible_color in unique_cylinders:
        if possible_color in aa:
            return possible_color
    return 'others'
def transform_owner_count(aa):
    if aa >= 5:
        return 'others'
    else:
        return aa
suv['interior_color'] = suv['interior_color'].apply(lambda x: transform_color(x))
suv['engine_cylinders'] = suv['engine_cylinders'].apply(lambda x: transform_engine_cylinders(x))
suv['owner_count'] = suv['owner_count'].apply(lambda x: transform_owner_count(x))
unique_colors = [color.upper() for color in suv['listing_color'].unique()]
unique_interior_color = [color.upper() for color in suv['interior_color'].unique()] 
unique_cylinders = suv['engine_cylinders'].unique()
unique_fuel_type = suv['fuel_type'].unique()
unique_transmission = suv['transmission'].unique()
unique_wheel_system = suv['wheel_system'].unique()
unique_make_name = suv['make_name'].unique()
unique_owner_count = suv['owner_count'].unique()
unique_maximum_seating = suv['maximum_seating'].unique()
categorical_cols = suv.columns[(17-1):(28-1+1)]
suv_dummies = pd.get_dummies(suv, columns=categorical_cols)
suv_dum_backup = suv_dummies.copy()
suv_dummies.drop('engine_cylinders_H4', axis=1, inplace=True)
suv_dummies.rename(columns={'make_name_Subaru': 'make_name_Subaru (engine_cylinders_H4)'}, inplace=True)
suv_dum_backup.drop('engine_cylinders_H4', axis=1, inplace=True)
suv_dum_backup.rename(columns={'make_name_Subaru': 'make_name_Subaru (engine_cylinders_H4)'}, inplace=True)
with open('index.json', 'r') as f:
    loaded_columns = json.load(f) 
columns_keep =  pd.Index(loaded_columns)
suv_dummies = suv_dummies[columns_keep[2:102]]
exterior_color_dummies = [col for col in suv_dummies.columns if 'listing_color' in col]
interior_color_dummies = [col for col in suv_dummies.columns if 'interior_color' in col]
fuel_type_dummies = [col for col in suv_dummies.columns if 'fuel_type' in col]
engine_cylinders_dummies = [col for col in suv_dummies.columns if 'engine_cylinders' in col]
transmission_dummies = [col for col in suv_dummies.columns if 'transmission' in col]
wheel_system_dummies = [col for col in suv_dummies.columns if 'wheel_system' in col]
make_name_dummies = [col for col in suv_dummies.columns if 'make_name' in col]
owner_count_dummies = [col for col in suv_dummies.columns if 'owner_count' in col]
maximum_seating_dummies = [col for col in suv_dummies.columns if 'maximum_seating' in col]
franchise_dealer_dummies = [col for col in suv_dummies.columns if 'franchise_dealer' in col]
isCab_dummies = [col for col in suv_dummies.columns if 'isCab' in col]
has_accidents_dummies = [col for col in suv_dummies.columns if 'has_accidents' in col]
for fuel_type_col in exterior_color_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['interior_color_BLACK']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['interior_color_BLACK'] ^ temp)/sum(suv_dum_backup['interior_color_BLACK'])
    k1 = 0.000001/len(interior_color_dummies)
    k2 = 0.000001/len(exterior_color_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in interior_color_dummies:
        temp = suv_dum_backup['listing_color_BLACK'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['listing_color_BLACK']^temp))/sum(suv_dum_backup['listing_color_BLACK'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp           
for fuel_type_col in fuel_type_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['wheel_system_AWD']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['wheel_system_AWD'] ^ temp)/sum(suv_dum_backup['wheel_system_AWD'])
    k1 = 0.000001/len(wheel_system_dummies)
    k2 = 0.000001/len(fuel_type_dummies)
    for engine_cyl_col in wheel_system_dummies:
        temp = suv_dum_backup['fuel_type_Gasoline'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['fuel_type_Gasoline']^temp))/sum(suv_dum_backup['fuel_type_Gasoline'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp         
for fuel_type_col in has_accidents_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['engine_cylinders_I4']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['engine_cylinders_I4'] ^ temp)/sum(suv_dum_backup['engine_cylinders_I4'])
    k1 = 0.000001/len(engine_cylinders_dummies)
    k2 = 0.000001/len(has_accidents_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in engine_cylinders_dummies:
        temp = suv_dum_backup['has_accidents_False'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['has_accidents_False']^temp))/sum(suv_dum_backup['has_accidents_False'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp
for fuel_type_col in has_accidents_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['wheel_system_AWD']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['wheel_system_AWD'] ^ temp)/sum(suv_dum_backup['wheel_system_AWD'])
    k1 = 0.000001/len(wheel_system_dummies)
    k2 = 0.000001/len(has_accidents_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in wheel_system_dummies:
        temp = suv_dum_backup['has_accidents_False'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['has_accidents_False']^temp))/sum(suv_dum_backup['has_accidents_False'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp         
for fuel_type_col in has_accidents_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['transmission_A']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['transmission_A'] ^ temp)/sum(suv_dum_backup['transmission_A'])
    k1 = 0.000001/len(transmission_dummies)
    k2 = 0.000001/len(has_accidents_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in transmission_dummies:
        temp = suv_dum_backup['has_accidents_False'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['has_accidents_False']^temp))/sum(suv_dum_backup['has_accidents_False'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp       
for fuel_type_col in has_accidents_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['fuel_type_Gasoline']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['fuel_type_Gasoline'] ^ temp)/sum(suv_dum_backup['fuel_type_Gasoline'])
    k1 = 0.000001/len(fuel_type_dummies)
    k2 = 0.000001/len(has_accidents_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in fuel_type_dummies:
        temp = suv_dum_backup['has_accidents_False'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['has_accidents_False']^temp))/sum(suv_dum_backup['has_accidents_False'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp 
for fuel_type_col in has_accidents_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['make_name_Jeep']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['make_name_Jeep'] ^ temp)/sum(suv_dum_backup['make_name_Jeep'])
    k1 = 0.000001/len(make_name_dummies)
    k2 = 0.000001/len(has_accidents_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in make_name_dummies:
        temp = suv_dum_backup['has_accidents_False'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['has_accidents_False']^temp))/sum(suv_dum_backup['has_accidents_False'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp       
for fuel_type_col in has_accidents_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['owner_count_1.0']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['owner_count_1.0'] ^ temp)/sum(suv_dum_backup['owner_count_1.0'])
    k1 = 0.000001/len(owner_count_dummies)
    k2 = 0.000001/len(has_accidents_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in owner_count_dummies:
        temp = suv_dum_backup['has_accidents_False'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['has_accidents_False']^temp))/sum(suv_dum_backup['has_accidents_False'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp        
for fuel_type_col in isCab_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['engine_cylinders_I4']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['engine_cylinders_I4'] ^ temp)/sum(suv_dum_backup['engine_cylinders_I4'])
    k1 = 0.000001/len(engine_cylinders_dummies)
    k2 = 0.000001/len(isCab_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in engine_cylinders_dummies:
        temp = suv_dum_backup['isCab_False'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['isCab_False']^temp))/sum(suv_dum_backup['isCab_False'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp 
for fuel_type_col in isCab_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['wheel_system_AWD']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['wheel_system_AWD'] ^ temp)/sum(suv_dum_backup['wheel_system_AWD'])
    k1 = 0.000001/len(wheel_system_dummies)
    k2 = 0.000001/len(isCab_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in wheel_system_dummies:
        temp = suv_dum_backup['isCab_False'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['isCab_False']^temp))/sum(suv_dum_backup['isCab_False'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp          
for fuel_type_col in isCab_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['transmission_A']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['transmission_A'] ^ temp)/sum(suv_dum_backup['transmission_A'])
    k1 = 0.000001/len(transmission_dummies)
    k2 = 0.000001/len(isCab_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in transmission_dummies:
        temp = suv_dum_backup['isCab_False'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['isCab_False']^temp))/sum(suv_dum_backup['isCab_False'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp        
for fuel_type_col in isCab_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['fuel_type_Gasoline']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['fuel_type_Gasoline'] ^ temp)/sum(suv_dum_backup['fuel_type_Gasoline'])
    k1 = 0.000001/len(fuel_type_dummies)
    k2 = 0.000001/len(isCab_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in fuel_type_dummies:
        temp = suv_dum_backup['isCab_False'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['isCab_False']^temp))/sum(suv_dum_backup['isCab_False'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp 
for fuel_type_col in isCab_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['make_name_Jeep']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['make_name_Jeep'] ^ temp)/sum(suv_dum_backup['make_name_Jeep'])
    k1 = 0.000001/len(make_name_dummies)
    k2 = 0.000001/len(isCab_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in make_name_dummies:
        temp = suv_dum_backup['isCab_False'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['isCab_False']^temp))/sum(suv_dum_backup['isCab_False'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp 
for fuel_type_col in isCab_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['owner_count_1.0']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['owner_count_1.0'] ^ temp)/sum(suv_dum_backup['owner_count_1.0'])
    k1 = 0.000001/len(owner_count_dummies)
    k2 = 0.000001/len(isCab_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in owner_count_dummies:
        temp = suv_dum_backup['isCab_False'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['isCab_False']^temp))/sum(suv_dum_backup['isCab_False'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp 
    for fuel_type_col in make_name_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['owner_count_1.0']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['owner_count_1.0'] ^ temp)/sum(suv_dum_backup['owner_count_1.0'])
    k1 = 0.000001/len(owner_count_dummies)
    k2 = 0.000001/len(make_name_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in owner_count_dummies:
        temp = suv_dum_backup['make_name_Jeep'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['make_name_Jeep']^temp))/sum(suv_dum_backup['make_name_Jeep'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp  
for fuel_type_col in exterior_color_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['wheel_system_AWD']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['wheel_system_AWD'] ^ temp)/sum(suv_dum_backup['wheel_system_AWD'])
    k1 = 0.1/len(wheel_system_dummies)
    k2 = 0.1/len(exterior_color_dummies)
    if c1<=k1 or c1>=(1-0.1) or c2<=k2 or c2>=(1-0.1):
        continue 
    for engine_cyl_col in wheel_system_dummies:
        temp = suv_dum_backup['listing_color_BLACK'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['listing_color_BLACK']^temp))/sum(suv_dum_backup['listing_color_BLACK'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.1) or c2<=k2 or c2>=(1-0.1):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.1) or c2<=k2 or c2>=(1-0.1):
            continue 
        suv_dummies[interaction_col_name] =  temp  
for fuel_type_col in exterior_color_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['engine_cylinders_I4']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['engine_cylinders_I4'] ^ temp)/sum(suv_dum_backup['engine_cylinders_I4'])
    k1 = 0.000001/len(engine_cylinders_dummies)
    k2 = 0.000001/len(exterior_color_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in engine_cylinders_dummies:
        temp = suv_dum_backup['listing_color_BLACK'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['listing_color_BLACK']^temp))/sum(suv_dum_backup['listing_color_BLACK'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp  
for fuel_type_col in exterior_color_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['fuel_type_Gasoline']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['fuel_type_Gasoline'] ^ temp)/sum(suv_dum_backup['fuel_type_Gasoline'])
    k1 = 0.000001/len(fuel_type_dummies)
    k2 = 0.000001/len(exterior_color_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in fuel_type_dummies:
        temp = suv_dum_backup['listing_color_BLACK'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['listing_color_BLACK']^temp))/sum(suv_dum_backup['listing_color_BLACK'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp  
for fuel_type_col in exterior_color_dummies:
    temp = suv_dummies[fuel_type_col] * suv_dum_backup['transmission_A']
    c1 = sum(suv_dum_backup[fuel_type_col]^temp)/sum(suv_dum_backup[fuel_type_col])
    c2 = sum(suv_dum_backup['transmission_A'] ^ temp)/sum(suv_dum_backup['transmission_A'])
    k1 = 0.000001/len(transmission_dummies)
    k2 = 0.000001/len(exterior_color_dummies)
    if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
        continue 
    for engine_cyl_col in transmission_dummies:
        temp = suv_dum_backup['listing_color_BLACK'] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup['listing_color_BLACK']^temp))/sum(suv_dum_backup['listing_color_BLACK'])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        interaction_col_name = f'{fuel_type_col}_X_{engine_cyl_col}'
        temp = suv_dummies[fuel_type_col] * suv_dummies[engine_cyl_col]
        c1 = (sum(suv_dum_backup[fuel_type_col]^temp))/sum(suv_dum_backup[fuel_type_col])
        c2 = (sum(suv_dum_backup[engine_cyl_col]^temp))/sum(suv_dum_backup[engine_cyl_col])
        if c1<=k1 or c1>=(1-0.000001) or c2<=k2 or c2>=(1-0.000001):
            continue 
        suv_dummies[interaction_col_name] =  temp  
suv_dummies1 = suv_dummies[columns_keep[2:]].copy()
df_sample = pd.DataFrame(suv_dummies1)
for column in df_sample.columns[1:15]:
    df_sample[column] = df_sample[column].astype(str).str.extract(r'(\d+\.?\d*)').astype(float)
df_sample.to_csv('suvtrain_reduced.csv', index=False)