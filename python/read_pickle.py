import pickle
#dd = {}

with open('db.pickle', 'rb') as f:
    dd = pickle.load(f)
    for key, value in dd.items():
        print('Observation: %r \nData: %r \n' %(key, value))

