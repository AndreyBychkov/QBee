class EvaluationStatistics:
    def __init__(self, steps: int, depth: int, method_name: str):
        self.steps = steps
        self.method_name = method_name
        self.depth = depth

    def __repr__(self):
        return '\n'.join([
            f"steps: {self.steps}",
            f"Method's name: {self.method_name}",
            f"depth: {self.depth}"
        ])
