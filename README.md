## Introduction
Possibility of a decentralized controller and associated communication savings is investigated for speed control of a string of moving vehicles using sparse Linear Quadratic Regulator (LQR).

## Abstract
With self-driving vehicles on the horizon, the control of vehicular formations is becoming increasingly lucrative. Of particular interest is a 1D formation wherein the vehicles move with a constant speed at fixed relative distances. For a small number of vehicles in the formation, all vehicles can easily communicate with one another and the problem is relatively easy to solve. However, as the number of vehicles in the formation increases, communications between individual vehicles becomes more difficult. A matlab routine was developed to solve the sparse LQR problem based on Alternating Direction Method of Multipliers. It was observed that the sparsity pattern of feedback matrix, F was near diagonal with less than 50 % non-zero elements. This meant that most vehicles will communicate with immediate or near immediate neighbours and still achieve a comparable performance.

*More details can be found in "Project_Report.pdf" file.*

*Current students: Please avoid plagiarizing this code as it results in violation of NC State's academic code of conduct.*
