
from sample_metadata.apis import ParticipantApi,SampleApi
from collections import Counter

participants = ParticipantApi().get_participants('tob-wgs')
 
print(participants[1])

age = []
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

master = [autoimmune_disease]
for i in participants: 
    try: 
        age.append(int(i['meta']['age']))
    except: 
        age.append(None)
    try: 
        autoimmune_disease.append((i['meta']['autoimmune_disease'].strip()))
    except: 
        autoimmune_disease.append(None)
    try: 
        diabetes_type1.append(i['meta']['diabetes_type1'].strip())
    except: 
        diabetes_type1.append(None)
    try: 
        diabetes_type2.append(i['meta']['diabetes_type2'].strip())
    except:
        diabetes_type2.append(None)
    try:
        rheumatoid_arthritis.append(i['meta']['rheumatoid_arthritis'].strip())      
    except: 
        rheumatoid_arthritis.append(None)
    try: 
        ulcerativecolitis.append(i['meta']['ulcerativecolitis'].strip())
    except: 
        ulcerativecolitis.append(None)
    try: 
        hyperthyroidism.append(i['meta']['hyperthyroidism'].strip())
    except: 
        hyperthyroidism.append(None)
    try: 
        hypertension.append(i['meta']['hypertension'].strip())
    except: 
        hypertension.append(None)
    try: 
        hypercholesterolaemia.append(i['meta']['hypercholesterolaemia'].strip())
    except: 
        hypercholesterolaemia.append(None)
    try: 
        cancer.append(i['meta']['cancer'].strip())
    except: 
        cancer.append(None)
    try: 
        eye_disease.append(i['meta']['eye_disease'].strip())
    except: 
        eye_disease.append(None)
    try: 
        osteoporosis.append(i['meta']['osteoporosis'].strip())
    except: 
        osteoporosis.append(None)
    try: 
        copd.append(i['meta']['copd'].strip())
    except: 
        copd.append(None)
    try: 
        statin.append(i['meta']['statin'].strip())
    except: 
        statin.append(None)
    try: 
        ace_inhibitor.append(i['meta']['ace_inhibitor'].strip())
    except: 
        ace_inhibitor.append(None)
    try: 
        angiotensinreptorblocker.append(i['meta']['angiotensinreptorblocker'].strip())
    except: 
        angiotensinreptorblocker.append(None)
    try: 
        calciumchannelblocker.append(i['meta']['calciumchannelblocker'].strip())
    except: 
        calciumchannelblocker.append(None)
    try: 
        betablocker.append(i['meta']['betablocker'].strip())
    except: 
        betablocker.append(None)
    try: 
        betaagonist.append(i['meta']['betaagonist'].strip())
    except: 
        betaagonist.append(None)
    try: 
        diuretic.append(i['meta']['diuretic'].strip())
    except: 
        diuretic.append(None)
    try: 
        oral_hypoglycaemic.append(i['meta']['oral_hypoglycaemic'].strip())
    except: 
        oral_hypoglycaemic.append(None)
    try: 
        insulin.append(i['meta']['insulin'].strip())
    except: 
        insulin.append(None)
    try: 
        ocp.append(i['meta']['ocp'].strip())
    except: 
        ocp.append(None)
    try: 
        paracetamol.append(i['meta']['paracetamol'].strip())
    except: 
        paracetamol.append(None)
    try: 
        aspirin.append(i['meta']['aspirin'].strip())
    except: 
        aspirin.append(None)
    try: 
        colchicine.append(i['meta']['colchicine'].strip())
    except: 
        colchicine.append(None)
    try: 
        ppi.append(i['meta']['ppi'].strip())
    except: 
        ppi.append(None)
    try: 
        thyroxine.append(i['meta']['thyroxine'].strip())
    except: 
        thyroxine.append(None)
    try: 
        smoking_status.append(i['meta']['smoking_status'].strip())
    except: 
        smoking_status.append(None)

master  = [autoimmune_disease,diabetes_type1,diabetes_type2,rheumatoid_arthritis,ulcerativecolitis,hyperthyroidism,hypertension,hypercholesterolaemia,cancer,
           eye_disease,osteoporosis,copd,statin,ace_inhibitor,angiotensinreptorblocker,calciumchannelblocker,betablocker,betaagonist,diuretic,oral_hypoglycaemic,insulin,
           ocp,paracetamol,aspirin,colchicine,ppi,thyroxine,smoking_status]
for m in master: 
    m_string= [i for i, j in locals().items() if j == m][0]
    print(f'{m_string}: {Counter(m)}')
    