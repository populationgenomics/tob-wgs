
from sample_metadata.apis import ParticipantApi,SampleApi

participants = ParticipantApi().get_participants('tob-wgs')
 
print(participants[1])

age = [9]
autoimmune_disease = []
diabetes_type1 =[]
diabetes_type2=[]
rheumatoid_arthritis= []
ulcerativecolitis =[]
hyperthyroidism=[]
autoimmune_disease_other=[]
hypertension=[]
hypercholesterolaemia=[]
cancer=[]
eye_disease=[]
osteoporosis=[]
copd=[]
statin=[]
ace_inhibitor=[]
angiotensinreptorblocker=[]
calciumchannelblocker=[]
betablocker=[]
betaagonist=[]
diuretic=[]
oral_hypoglycaemic=[]
insulin=[]
ocp=[]
paracetamol=[]
aspirin=[]
colchicine=[]
ppi=[]
thyroxine=[]
smoking_status=[]

for i in participants: 
    try: 
        age.append(int(i['meta']['age']))
    except: 
        print(i['external_id'])
        age.append(None)
   


