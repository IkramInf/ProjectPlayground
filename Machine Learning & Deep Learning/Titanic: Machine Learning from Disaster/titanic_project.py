#!/usr/bin/env python
# coding: utf-8

# In[73]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import accuracy_score


# In[74]:


train = pd.read_csv('/home/ikram2718/Desktop/Python_program/datasets/titanic/train.csv')
test = pd.read_csv('/home/ikram2718/Desktop/Python_program/datasets/titanic/test.csv')
#y_test = pd.read_csv('/home/ikram2718/Desktop/Python_program/datasets/titanic/titanic/gender_submission.csv')


# In[75]:


train.columns


# In[76]:


#y_test = y_test['Survived']


# In[77]:


def count_param(feature):
    survived = train[train['Survived']==1][feature].value_counts()
    dead = train[train['Survived']==0][feature].value_counts()
    print(survived, dead)


# In[78]:


train_data = train.drop(['PassengerId','Name','Ticket','Cabin'],axis=1)
test_data = test.drop(['PassengerId','Name','Ticket','Cabin'],axis=1)


# In[79]:


train_data.info()


# In[80]:


test_data.info()


# In[81]:


train_data['Embarked'] = train_data['Embarked'].fillna(value='S')


# In[111]:


train_data['Age'].fillna(value=train_data.groupby(['Parch','SibSp','Pclass'])['Age'].transform('mean'), inplace=True)
test_data['Age'].fillna(value=test_data.groupby(['Parch','SibSp'])['Age'].transform('mean'), inplace=True)


# In[112]:


train_data = train_data.interpolate()
test_data = test_data.interpolate()


# In[113]:


train_data.info()


# In[114]:


test_data.info()


# In[115]:


from sklearn.preprocessing import LabelEncoder
label_encoder = LabelEncoder()
#sex = ['male','female']
train_data['Sex'] = label_encoder.fit_transform(train_data['Sex'])
test_data['Sex'] = label_encoder.fit_transform(test_data['Sex'])


# In[116]:


#embarked = ['S','Q','C']
train_data['Embarked'] = label_encoder.fit_transform(train_data['Embarked'])
test_data['Embarked'] = label_encoder.fit_transform(test_data['Embarked'])


# In[117]:


X = train_data.drop(['Survived'],axis=1)
y = train_data['Survived']


# In[118]:


print(type(X))


# In[119]:


from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.20, random_state = 1)


# In[120]:


from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
k_fold = KFold(n_splits=10, shuffle=True, random_state=1)


# In[121]:


from sklearn.tree import DecisionTreeClassifier
model = DecisionTreeClassifier(criterion='entropy', max_depth=3, min_samples_leaf=2, random_state=1)
score = cross_val_score(model, X_train, y_train, cv=k_fold, n_jobs=1, scoring='accuracy')
print(score)
print(np.mean(score))


# In[122]:


from sklearn.neighbors import KNeighborsClassifier
clf = KNeighborsClassifier(n_neighbors = 2)
score = cross_val_score(clf, X_train, y_train, cv=k_fold, n_jobs=1, scoring='accuracy')
print(score)
print(np.mean(score))


# In[123]:


from sklearn.ensemble import RandomForestClassifier
clf_r = RandomForestClassifier(n_estimators=31, criterion='entropy', max_depth=10)
score = cross_val_score(clf_r, X_train, y_train, cv=k_fold, n_jobs=1, scoring='accuracy')
print(score)
round(np.mean(score)*100, 2)


# In[124]:


model.fit(X,y)
y_pred = model.predict(test_data)


# In[125]:


submission = pd.DataFrame({"PassengerId" : test['PassengerId'], 'Survived': y_pred})
submission.to_csv('Titanic_project_8th_submission.csv', index=False)


# In[126]:


clf_r.fit(X_train,y_train)
y_pred = clf_r.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)
print(accuracy)


# In[127]:


model.fit(X_train,y_train)
y_pred = model.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)
print(accuracy)


# In[110]:


clf.fit(X_train, y_train)
y_p = clf.predict(X_test)


# In[100]:


acr = accuracy_score(y_test, y_p)
print(acr)
