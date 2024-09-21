# Hand-Shape-Approximation-with-Clamped-Cubic-Spline
### Purpose/Objective
To approximate the shape of my left hand using the Clamped Cubic Spline Algorithm and numerical differentiation techniques.

### Methodology
Collected key hand outline points and divided each finger into two curves (10 curves total). Implemented 3-point forward/backward differentiation to estimate slopes at endpoints and used the Clamped Cubic Spline Algorithm to create smooth curves. Adjusted the data for improved accuracy.

### Tools
- MATLAB: For implementing the algorithms and plotting the results.
- 3-Point Differentiation: For calculating endpoint slopes.
- Clamped Cubic Spline: For creating smooth curves.

### Results
- The spline curves closely matched the hand’s shape, with some minor differences, especially in the middle finger, due to data collection limitations.
- The Clamped Cubic Spline Algorithm and 3-point differentiation successfully approximated the hand’s shape, producing a smooth and accurate plot.
