-- t1: table name of your data group A, use NTILE(4) as an example 
-- t3: temp data table with calculated y_mean, x thresholds, save the results as t3
-- t4: table name of your data group B

WITH t2 AS (SELECT t1_id, x, NTILE(4) OVER (ORDER BY x) 
   AS quantile FROM t1)
   SELECT quantile,
   AVG(t2.y) as y_mean,
   MAX(t2.x) as x_max,
   MIN(t2.x) as x_min
   FROM t2 
   GROUP BY quantile;

-- based on the result from above say table3 (t3), you know the x_max1-4
-- assign the quntile group for table 4 (your group B data)

SELECT x, y
(CASE WHEN x <= x_max1 THEN 1
   WHEN x > x_max1 and x <= x_max2 THEN 2
   WHEN x > x_max2 and x <= x_max3 THEN 3
   WHEN x > x_max3 and x <= x_max4 THEN 4
   ELSE 5 END) as quantile
from t4;

select t4.y, t4.quantile, t3.y_mean,
(case when y>y_mean THEN 1 ELSE 0 END) as select_label
left join t3 on t3.quantile=t4.quantile; 


-- #########################################################
WITH t5 AS (SELECT x, y
(CASE WHEN x <= x_max1 THEN 1
   WHEN x > x_max1 and x <= x_max2 THEN 2
   WHEN x > x_max2 and x <= x_max3 THEN 3
   WHEN x > x_max3 and x <= x_max4 THEN 4
   ELSE 5 END) AS quantile
FROM t4) 
SELECT t5.y, t5.quantile, t3.y_mean,
(case when y>y_mean THEN 1 ELSE 0 END) as select_label from t3,t5
where t3.quantile=t5.quantile; 
