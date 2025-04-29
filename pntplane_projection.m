%input(nx3), parameter=[A;B;C;D]  即平面方程为Ax+By+Cz+D=0

function [project_pnt] = pntplane_projection(input,parameter)    %该方法效果最好
    A=parameter(1);B=parameter(2);C=parameter(3);D=parameter(4);
    a11=B^2+C^2;a12=-A*B;  a13=-A*C;
    a21=-A*B;  a22=A^2+C^2;a23=-B*C;
    a31=-A*C;  a32=-B*C;  a33=A^2+B^2;
    b1=-A*D;    b2=-B*D;   b3=-C*D;
    a_matrix=[a11,a12,a13;a21,a22,a23;a31,a32,a33];
    b_matrix=[b1;b2;b3];
    project_pnt=a_matrix*input'+b_matrix;
    
