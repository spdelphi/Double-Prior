Our algorithm: adsu_pq.m
Easy example to run our algorithm: demo_easy.m
Full comparision in DC1 of our paper: demo_DC1.m
Demo for real dataset: demo_real_datasets.m
SVM Classifier( one way to produce prior): classification.py
By the classification result, we can use the following code to calculate prior:
w = ones(bands,nCol*nRow);
w(sub2ind(size(w), classification_result, 1:length(classification_result))) = 0;
w = normWeight(w);