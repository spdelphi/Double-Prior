# Double prior model

Our algorithm: adsu_pq.m

Easy example to run our algorithm: demo_easy.m

Full comparision in DC1 of our paper: demo_DC1.m

Demo for real dataset: demo_real_datasets.m

SVM Classifier( one way to generate prior): classification.py

With the classification result, we can use the following code to calculate the prior:

w = ones(bands,nCol*nRow);

w(sub2ind(size(w), classification_result, 1:length(classification_result))) = 0;

w = normWeight(w);



If you use the codes, please cite the following paper. Thank you very much.

R. Wu, W. Luo, L. Gao, L. Ren and H. Gao, "Endmember Selection With Adaptive Double Prior Model," in IEEE Transactions on Geoscience and Remote Sensing, vol. 64, pp. 1-22, 2026, Art no. 5508122, doi: 10.1109/TGRS.2026.3664867.
