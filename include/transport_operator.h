#ifndef TRANSPORT_OPERATOR_H_
#define TRANSPORT_OPERATOR_H_

class TransportOperator {
public:
	virtual void CreateTallies() = 0;
	virtual void SimulateBatches(const int& num_batches) = 0;
};

#endif