function data = cutextreme(data,tail)
lower = quantile(data, tail);
upper = quantile(data, 1-tail);
data = data(data >= lower & data <= upper);