## Changing the traversal order instead of transforming the pixmap
Row-by-row traversal ("mirroring" -v, -h, and -vh (same as r180)) is free (for speed) and cache-friendly.

The cost of traversing column-by-column ("rotating" 90, 90v, 90h, -90) is implementation dependent.
