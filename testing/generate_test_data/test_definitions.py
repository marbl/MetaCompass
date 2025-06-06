from shuffle import shuffle
from substring import substring


class tests:

    available_test_definitions = {
        shuffle.name: shuffle.apply,
        substring.name: substring.apply,
    }
    available_test_names = list(available_test_definitions.keys())

    def apply(self):
        print("error: apply not initialized")
        return "error: apply not set"

    def get_apply_function(self):
        if not self.valid:
            return

        self.apply = self.available_test_definitions[self.test_name]
        return

    def is_valid(self):
        if self.test_name not in self.available_test_names:
            return False

        return True

    def __init__(self, args):
        self.test_name = args.test_name
        self.valid = self.is_valid()
        self.seed = args.seed
        if not self.valid:
            return

        print("initializing test type " + str(self.test_name))

        self.get_apply_function()
        return
