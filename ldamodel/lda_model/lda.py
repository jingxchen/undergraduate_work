# -*- coding: utf-8 -*-

import sys
import os
import pandas as pd
import numpy as np
from gensim import corpora, models, similarities
import logging
from collections import defaultdict
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize
from nltk.corpus import brown
from nltk.stem.lancaster import LancasterStemmer
import pyLDAvis.gensim
from matplotlib import pyplot as plt
st = LancasterStemmer()


#打开文章
path1 = 'C:/Users/95260/Desktop/'
path2 = 'C:/Users/95260/Desktop'
path3 = 'C:/Users/95260/Desktop'
#path1 = sys.argv[1]
#path2 = sys.argv[2]
#path3 = sys.argv[3]
#path1 = 论文文件的存放路径
#path2 = 输出图片的存放路径
#path3 = 论文年份的存放路径


files= os.listdir(path1)
article=[]
for file in files:
     if not os.path.isdir(file):
          f = open(path1+"/"+file, 'r', encoding = 'UTF-8')
          line=f.readline()
          while(line):
              article.append(line)
              line=f.readline()
print(len(article))
f.close()
texts_tokenized = [[word.lower() for word in word_tokenize(document)] for document in article]
#去停用词
stoplist = stopwords.words('english')
useful_words = [[word for word in document if not word in stoplist] for document in texts_tokenized]
#去标点
english_punctuations = [',', '.', ':', ';', '?', '(', ')', '[', ']', '&', '!', '*', '@', '#', '$', '%']
useful_words_nopunc = [[word for word in document if not word in english_punctuations] for document in useful_words]
#词干化
#stemmer = LancasterStemmer()
#texts = [[stemmer.stem(word) for word in docment] for docment in useful_words_nopunc]
texts = useful_words_nopunc
dictionary = corpora.Dictionary(texts)
corpus = [dictionary.doc2bow(text) for text in texts]
tfidf = models.TfidfModel(corpus)
corpus_tfidf = tfidf[corpus]

#训练LDA模型
lda = models.LdaModel(corpus_tfidf, id2word=dictionary, num_topics=5,passes = 10)

#可视化

data1=pyLDAvis.gensim.prepare(lda,corpus,dictionary)
#pyLDAvis.show(data1,open_browser=True)
pyLDAvis.save_html(data1,path2+'term freq.html')
pyLDAvis.save_json(data1,path2+'term freq.json')

#from sklearn.externals import joblib
#joblib.dump(lda,'D:/kg/111/model_save.m')
#lda=joblib.load('D:/kg/111/model_save.m')

def takeSecond(elem):
    return elem[1]

topic1=[1]*len(article)
lgt = lda.get_document_topics(bow=corpus_tfidf)
for i in range(len(article)):
    lgti=lgt[i]
    lgti.sort(key=takeSecond,reverse=True)
    topic1[i]=lgti[0][0]

f=open(path3+'year.txt')
numbers = f.readlines()
array = []
for i in numbers:
    i = i.replace("\n","")
    i = i.split(" ")
    i = int(i[0])
    array.append(i)

years = array[0:len(article)]
count = np.zeros((max(years)-min(years)+1, 10))
for i in range(len(years)):
    count[years[i] - min(years)][topic1[i]] = count[years[i] - min(years)][topic1[i]] + 1

color_sequence = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c',
                  '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5',
                  '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f',
                  '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']

x =  range(min(years),max(years)+1)
plt.style.use('fivethirtyeight')
for i in range(10):
    plt.title('Trend analysis for topic '+str(i))
    plt.plot(x, count[:,i],label="topic"+str(i),color=color_sequence[i],marker='o',mec=color_sequence[i+2],mfc=color_sequence[i+3])
    plt.text(2015.5, count[-1,i], "topic " + str(i), fontsize=14, color='black')
    plt.xlabel('Year')
    plt.ylabel('Numbers of documents')
 #   plt.xlim(2004, 2017)
  #  plt.ylim(0,12000)
    plt.savefig(path2+'topic '+str(i)+'.png', bbox_inches='tight')
    plt.show()


for i in range(10):
    plt.plot(x, count[:,i],label="topic"+str(i),color=color_sequence[i],marker='o',mec=color_sequence[i+2],mfc=color_sequence[i+3])
    plt.text(2015.5, count[-1,i], "topic " + str(i), fontsize=10, color='black')

plt.title('Trend analysis for all topics ')
plt.xlabel('Year')
plt.ylabel('Numbers of documents')
#plt.xlim(2004, 2017)
#plt.ylim(0,12000)
plt.legend(loc=2)
plt.savefig(path2+'all topics.png', bbox_inches='tight')
plt.show()

