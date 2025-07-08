-- Check where in the code this is used
-- It may be able to be deleted 

DELETE FROM pairs
    WHERE rowid NOT IN (
        SELECT MIN(rowid)
        FROM (
            SELECT rowid, 
                    LEAST(strain_1, strain_2) AS s1, 
                    GREATEST(strain_1, strain_2) AS s2
            FROM pairs
        )
        GROUP BY s1, s2
);
