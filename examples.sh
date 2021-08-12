# Run the algorithm
./build/manifold --input data/bathtub.obj --output results/bathtub_manifold.obj
./build/manifold --input data/bed.obj --output results/bed_manifold.obj
./build/manifold --input data/chair.off --output results/chair_manifold.obj
./build/manifold --input data/table.off --output results/table_manifold.obj
echo 'Check results in results folder. If you use meshlab, turn on double face rendering since faces from both sides of zero-volume structure could collide together.'
