# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 19:46:25 2018

@author: Administrator
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from nltk.tokenize import TweetTokenizer
import datetime
import lightgbm as lgb
from scipy import stats
from scipy.sparse import hstack, csr_matrix
from sklearn.model_selection import train_test_split, cross_val_score
from wordcloud import WordCloud
from collections import Counter
from nltk.corpus import stopwords
from nltk.util import ngrams
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.multiclass import OneVsRestClassifier
from sklearn.tree import DecisionTreeClassifier 
from sklearn.naive_bayes import BernoulliNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import RidgeClassifier
from sklearn.neighbors import NearestCentroid
pd.set_option('max_colwidth',400)

#data input
train = pd.read_csv('c:/Users/Administrator/Desktop/kw/train.csv',sep=',')
test = pd.read_csv('C:/Users/Administrator/Desktop/kw/test.csv',sep=',')
sub = pd.read_csv('C:/Users/Administrator/Desktop/kw/sampleSubmission.csv', sep=",")

#data analysis
train.head(10)
print('Average number of phrases per sentence in train is {0:.0f}.'.format(train.groupby('SentenceId')['Phrase'].count().mean()))
print('Average number of phrases per sentence in test is {0:.0f}.'.format(test.groupby('SentenceId')['Phrase'].count().mean()))

print('Number of phrases in train: {}. Number of sentences in train: {}.'.format(train.shape[0], len(train.SentenceId.unique())))
print('Number of phrases in test: {}. Number of sentences in test: {}.'.format(test.shape[0], len(test.SentenceId.unique())))

print('Average word length of phrases in train is {0:.0f}.'.format(np.mean(train['Phrase'].apply(lambda x: len(x.split())))))
print('Average word length of phrases in test is {0:.0f}.'.format(np.mean(test['Phrase'].apply(lambda x: len(x.split())))))

text4 = ' '.join(train.loc[train.Sentiment == 4, 'Phrase'].values)
text4_trigrams = [i for i in ngrams(text4.split(), 4)]
Counter(text4_trigrams).most_common(10)

text0 = ' '.join(train.loc[train.Sentiment == 0, 'Phrase'].values)
text0_trigrams = [i for i in ngrams(text0.split(), 4)]
Counter(text0_trigrams).most_common(10)

text4 = ' '.join(train.loc[train.Sentiment == 4, 'Phrase'].values)
text4_nostopword = [i for i in text4.split() if i not in stopwords.words('english')]
text4_nostopword_trigrams = [i for i in ngrams(text4_nostopword, 4)]
Counter(text4_nostopword_trigrams).most_common(10)

english_punctuations = [',', '.', ':', ';', '?', '(', ')', '[', ']', '&', '!', '*', '@', '#', '$', '%']
text4_nopunctuation = [i for i in text4.split() if i not in english_punctuations]
text4_nopunctuation_trigrams = [i for i in ngrams(text4_nopunctuation, 4)]
Counter(text4_nopunctuation_trigrams).most_common(10)

english_punctuations = [',', '.', ':', ';', '?', '(', ')', '[', ']', '&', '!', '*', '@', '#', '$', '%']
text0_nopunctuation = [i for i in text0.split() if i not in english_punctuations]
text0_nopunctuation_trigrams = [i for i in ngrams(text0_nopunctuation, 4)]
Counter(text0_nopunctuation_trigrams).most_common(10)
#seperate
tokenizer = TweetTokenizer()
#vectorizer = TfidfVectorizer(ngram_range=(1, 2), tokenizer=tokenizer.tokenize)
vectorizer = TfidfVectorizer(tokenizer=tokenizer.tokenize)
full_text = list(train['Phrase'].values) + list(test['Phrase'].values)
vectorizer.fit(full_text)
train_vectorized = vectorizer.transform(train['Phrase'])
test_vectorized = vectorizer.transform(test['Phrase'])
aaa=vectorizer.get_feature_names()
aaa.index('of')

sentiment = train['Sentiment']

#1
logreg = LogisticRegression()
ovr = OneVsRestClassifier(logreg)
of = ovr.fit(train_vectorized, sentiment);
of.score(train_vectorized,sentiment)
scores_o = cross_val_score(ovr, train_vectorized, sentiment, scoring='accuracy', n_jobs=-1, cv=3)
print('Cross-validation mean accuracy {0:.2f}%, std {1:.2f}.'.format(np.mean(scores_o) * 100, np.std(scores_o) * 100))
#Cross-validation mean accuracy 58.31%, std 0.07. fit.score: 0.67055619633474306

#2
svc = LinearSVC(dual=False)
sf = svc.fit(train_vectorized,sentiment);
sf.score(train_vectorized,sentiment)
scores_s = cross_val_score(svc, train_vectorized, sentiment, scoring='accuracy', n_jobs=-1, cv=3)
print('Cross-validation mean accuracy {0:.2f}%, std {1:.2f}.'.format(np.mean(scores_s) * 100, np.std(scores_s) * 100))
#Cross-validation mean accuracy 57.81%, std 0.54. fit.score: 0.72063949762911705

#3
dec = DecisionTreeClassifier()
df = dec.fit(train_vectorized, sentiment);
df.score(train_vectorized,sentiment)
scores_d = cross_val_score(dec, train_vectorized, sentiment, scoring='accuracy', n_jobs=-1, cv=3)
print('Cross-validation mean accuracy {0:.2f}%, std {1:.2f}.'.format(np.mean(scores_d) * 100, np.std(scores_d) * 100))
#Cross-validation mean accuracy 47.03%, std 3.38. fit.score: 0.99911572472126109

#4
ber = BernoulliNB()
bf = ber.fit(train_vectorized, sentiment);
bf.score(train_vectorized,sentiment)
scores_b = cross_val_score(ber, train_vectorized, sentiment, scoring='accuracy', n_jobs=-1, cv=3)
print('Cross-validation mean accuracy {0:.2f}%, std {1:.2f}.'.format(np.mean(scores_b) * 100, np.std(scores_b) * 100))
#Cross-validation mean accuracy 53.98%, std 0.62. fit.score: 0.62381135460720238

#5
ran = RandomForestClassifier()
rf = ran.fit(train_vectorized, sentiment);
rf.score(train_vectorized,sentiment)
scores_r = cross_val_score(ran, train_vectorized, sentiment, scoring='accuracy', n_jobs=-1, cv=3)
print('Cross-validation mean accuracy {0:.2f}%, std {1:.2f}.'.format(np.mean(scores_r) * 100, np.std(scores_r) * 100))
#Cross-validation mean accuracy 50.15%, std 2.98. fit.score: 0.96962706651287967
   
#7
rid = RidgeClassifier()
rf = rid.fit(train_vectorized, sentiment);
rf.score(train_vectorized,sentiment)
scores_r = cross_val_score(rid, train_vectorized, sentiment, scoring='accuracy', n_jobs=-1, cv=3)
print('Cross-validation mean accuracy {0:.2f}%, std {1:.2f}.'.format(np.mean(scores_r) * 100, np.std(scores_r) * 100))
#Cross-validation mean accuracy 57.91%, std 0.45 . fit.score:  0.69803280789439959
 
#8

#lstm
import numpy as np 
import pandas as pd 
import nltk
import os
import gc
from keras.preprocessing import sequence,text
from keras.preprocessing.text import Tokenizer
from keras.models import Sequential
from keras.layers import Dense,Dropout,Embedding,LSTM,Conv1D,GlobalMaxPooling1D,Flatten,MaxPooling1D,GRU,SpatialDropout1D,Bidirectional
from keras.callbacks import EarlyStopping
from keras.utils import to_categorical
from keras.losses import categorical_crossentropy
from keras.optimizers import Adam
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score,confusion_matrix,classification_report,f1_score
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")   

train = pd.read_csv('c:/Users/Administrator/Desktop/kw/train.csv',sep=',')    
print(train.shape)
train.head()
test = pd.read_csv('C:/Users/Administrator/Desktop/kw/test.csv',sep=',')
sub = pd.read_csv('C:/Users/Administrator/Desktop/kw/sampleSubmission.csv', sep=",")
test['Sentiment']=-999
df=pd.concat([train,test],ignore_index=True)
print(df.shape)
del train,test
gc.collect()

from nltk.tokenize import word_tokenize
from nltk import FreqDist
from nltk.stem import SnowballStemmer,WordNetLemmatizer
stemmer=SnowballStemmer('english')
lemma=WordNetLemmatizer()
from string import punctuation
import re

def clean_review(review_col):
    review_corpus=[]
    for i in range(0,len(review_col)):
        review=str(review_col[i])
        review=re.sub('[^a-zA-Z]',' ',review)
        #review=[stemmer.stem(w) for w in word_tokenize(str(review).lower())]
        review=[lemma.lemmatize(w) for w in word_tokenize(str(review).lower())]
        review=' '.join(review)
        review_corpus.append(review)
    return review_corpus
    
df['clean_review']=clean_review(df.Phrase.values)
df.head()    
df_train=df[df.Sentiment!=-999]
df_test=df[df.Sentiment==-999]
df_test.drop('Sentiment',axis=1,inplace=True)
del df
gc.collect()

train_text=df_train.clean_review.values
test_text=df_test.clean_review.values
target=df_train.Sentiment.values
y=to_categorical(target)
print(train_text.shape,target.shape,y.shape)

#cv
X_train_text,X_val_text,y_train,y_val=train_test_split(train_text,y,test_size=0.2,stratify=y,random_state=123)

all_words=' '.join(X_train_text)
all_words=word_tokenize(all_words)
dist=FreqDist(all_words)
num_unique_word=len(dist)

r_len=[]
for text in X_train_text:
    word=word_tokenize(text)
    l=len(word)
    r_len.append(l)
    
MAX_REVIEW_LEN=np.max(r_len)
MAX_REVIEW_LEN

max_features = num_unique_word
max_words = MAX_REVIEW_LEN
batch_size = 128
epochs = 3
num_classes=5

tokenizer = Tokenizer(num_words=max_features)
tokenizer.fit_on_texts(list(X_train_text))
X_train = tokenizer.texts_to_sequences(X_train_text)
X_val = tokenizer.texts_to_sequences(X_val_text)
X_test = tokenizer.texts_to_sequences(test_text)

#sequence padding?
X_train = sequence.pad_sequences(X_train, maxlen=max_words)
X_val = sequence.pad_sequences(X_val, maxlen=max_words)
X_test = sequence.pad_sequences(X_test, maxlen=max_words)

model1=Sequential()
model1.add(Embedding(max_features,100,mask_zero=True))
model1.add(LSTM(64,dropout=0.4, recurrent_dropout=0.4,return_sequences=True))
model1.add(LSTM(32,dropout=0.5, recurrent_dropout=0.5,return_sequences=False))
model1.add(Dense(num_classes,activation='softmax'))
model1.compile(loss='categorical_crossentropy',optimizer=Adam(lr=0.001),metrics=['accuracy'])
model1.summary()

history1=model1.fit(X_train, y_train, validation_data=(X_val, y_val),epochs=3, batch_size=batch_size, verbose=1)
#Cross-validation mean accuracy 67%