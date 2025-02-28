Repo for velia DB/Data Lake dashboards

Installation

git clone https://github.com/veliatx/orfdb.git
pip install -e ./orfdb

git clone https://github.com/veliatx/conservation.git
pip install -e ./conservation

git clone https://github.com/veliatx/dashboard.git
cd dashboard
pip install -r requirements.txt

Usage
#Add new orfs to the dashboard.  
Note that /efs must be mounted to access the signalP 6 weights

#Depends on ORF being present in ORFDB  
#Create input txt file with VTX IDs one id per line  
#dashboard_update --help for more info  
dashboard_update /path/to/vtx_input_file.txt /path/to/output/cache /path/to/data  

#Launches the streamlit app
streamlit run /path/to/sorf_app.py


