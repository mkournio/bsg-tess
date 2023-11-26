import logging

log = logging.getLogger("my-logger")
log.info("Hello, world")
ARGS = (3,4,5)

def test(x,*ARGS): 

	print x, ARGS
	log.warning(
                "You have passed both a grid of frequencies "
                "and min_frequency/maximum_frequency arguments; "
                "the latter will be ignored."
            )

print test(3, 3)


