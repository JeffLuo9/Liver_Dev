import pandas as pd
import numpy as np
import skimage as ski
import os
import glob
from pathlib import Path
os.chdir("F:/liver_image/")

def calculate_dapi_bright_analysis(dapi_raw_path, bright_path, dapi_path,sample_name ,output_file,cover_thre=0.5):
    import pandas as pd
    import numpy as np
    import skimage as ski

    dapi_raw_img = ski.io.imread(dapi_raw_path)
    bright_img = ski.io.imread(bright_path)
    dapi = pd.read_csv(dapi_path)

    y,x = np.where(bright_img!=0)
    bright = pd.DataFrame({'y':y,'x':x,'value':bright_img[y,x]})
    dapi['bright'] = bright_img[dapi['y'],dapi['x']]
    dapi['area'] = 1
    dapi['bright_area'] = dapi['bright'].copy()
    dapi_agg = dapi.groupby(['value']).agg({'bright': lambda x: x.value_counts().index[0],
                                        'bright_area': lambda x: x.value_counts().iloc[0],
                                        'area':'sum'
                                       }).reset_index()

    cover_thre = 0.5
    dapi_agg = dapi_agg[dapi_agg['bright']!=0].copy()
    dapi_agg['cover_percentage'] = dapi_agg['bright_area']/dapi_agg['area']
    dapi_agg = dapi_agg[dapi_agg['cover_percentage']>cover_thre].copy()
    dapi_dtb = dapi_agg

    tmp = np.zeros((max([bright['y'].max(),dapi_raw_img.shape[0]])+5,
                max([bright['x'].max(),dapi_raw_img.shape[1]])+5))
    tmp[0:dapi_raw_img.shape[0],0:dapi_raw_img.shape[1]] = dapi_raw_img
    bright['dapi_value'] = tmp[bright['y'],bright['x']]
    bright['area'] = 1
    bright['x_y'] = bright['x'].map(str) + '_' + bright['y'].map(str)
    dapi['x_y'] = dapi['x'].map(str) + '_' + dapi['y'].map(str)
    bright['dapi_area'] = bright['x_y'].isin(dapi['x_y'])
    bright_agg = bright.groupby(['value']).agg({'x':'mean',
                                            'y':'mean',
                                            'dapi_value':'mean',
                                            'dapi_area':'sum',
                                            'area':'sum'}).reset_index()
    bright_agg['dapi_area_percentage'] = bright_agg['dapi_area']/bright_agg['area']
    dapi_agg['count'] = 1
    dapi_agg_counts = dapi_agg.groupby(['bright']).agg({'count':'sum'}).reset_index()
    bright_agg['dapi_counts'] = bright_agg['value'].map(dict(zip(dapi_agg_counts['bright'],dapi_agg_counts['count'])))
    bright_agg = bright_agg.fillna(0)
    bright_agg['dapi_counts'] = bright_agg['dapi_counts'].astype(int)
    dapi_dtb.to_csv(output_file+sample_name+"_dtb.csv",index = False)
    bright_agg.to_csv(output_file+sample_name+"_nucle_info.csv",index=False)
    
    return 1

sample_list = ["Y00723K8", "Y00723J6", "Y00723E5", "Y00724A4", "Y00724A1","Y00723G8","Y00723H8","Y00723PA","Y00724A2","Y00788PC","Y00788N9"
               "A04224G4", "A04224G5", "A04224E6", "A04225F2", "Y00723H3","Y00724PD","Y00723M1"]

for sample_name in sample_list:
    bright_path = f"F:/liver_image/nac-{sample_name}_cp_masks.tif"
    dapi_raw_path = f"F:/liver_image/nac-dapi-{sample_name}.tif"
    dapi_path = f"F:/liver_image/nac-dapi-{sample_name}_new.csv"

    calculate_dapi_bright_analysis(
        dapi_raw_path=dapi_raw_path,
        bright_path=bright_path,
        dapi_path=dapi_path,
        sample_name=sample_name,
        output_file="F:/liver_image/dtb_files/"
    
    )
    print(sample_name)


