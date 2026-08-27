# EstimationGears Consumer Project

This disposable example verifies the installed CMake package and exported target.
Configure it after installing the main project:

```bash
cmake -S examples/template_consumer_project -B /tmp/estimation-gears-consumer \
  -DCMAKE_PREFIX_PATH=/path/to/estimation-gears-install
cmake --build /tmp/estimation-gears-consumer
/tmp/estimation-gears-consumer/estimation_gears_consumer
```
