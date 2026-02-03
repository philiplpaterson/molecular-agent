import { useState } from "react";
import { ToolCall } from "./MessageList";

interface ToolCallDisplayProps {
  toolCall: ToolCall;
}

const TOOL_ICONS: Record<string, string> = {
  validate_smiles: "✓",
  predict_properties: "📊",
  check_drug_likeness: "💊",
  similarity_search: "🔍",
  generate_molecules: "🧬",
};

const TOOL_LABELS: Record<string, string> = {
  validate_smiles: "Validate SMILES",
  predict_properties: "Predict Properties",
  check_drug_likeness: "Drug-Likeness Check",
  similarity_search: "Similarity Search",
  generate_molecules: "Generate Molecules",
};

export default function ToolCallDisplay({ toolCall }: ToolCallDisplayProps) {
  const [expanded, setExpanded] = useState(false);

  const icon = TOOL_ICONS[toolCall.name] || "🔧";
  const label = TOOL_LABELS[toolCall.name] || toolCall.name;

  const hasError = "error" in toolCall.result && toolCall.result.error != null;

  return (
    <div className="bg-white border border-gray-200 rounded-lg overflow-hidden text-sm">
      <button
        onClick={() => setExpanded(!expanded)}
        className="w-full px-3 py-2 flex items-center justify-between hover:bg-gray-50 transition-colors"
      >
        <div className="flex items-center gap-2">
          <span>{icon}</span>
          <span className="font-medium text-gray-700">{label}</span>
          {hasError && (
            <span className="text-xs bg-red-100 text-red-600 px-2 py-0.5 rounded">
              Error
            </span>
          )}
        </div>
        <svg
          className={`w-4 h-4 text-gray-400 transition-transform ${
            expanded ? "rotate-180" : ""
          }`}
          fill="none"
          stroke="currentColor"
          viewBox="0 0 24 24"
        >
          <path
            strokeLinecap="round"
            strokeLinejoin="round"
            strokeWidth={2}
            d="M19 9l-7 7-7-7"
          />
        </svg>
      </button>

      {expanded && (
        <div className="px-3 py-2 border-t border-gray-200 bg-gray-50 space-y-2">
          <div>
            <div className="text-xs font-medium text-gray-500 mb-1">Input</div>
            <pre className="text-xs bg-gray-100 p-2 rounded overflow-x-auto">
              {JSON.stringify(toolCall.arguments, null, 2)}
            </pre>
          </div>
          <div>
            <div className="text-xs font-medium text-gray-500 mb-1">Output</div>
            <pre
              className={`text-xs p-2 rounded overflow-x-auto ${
                hasError ? "bg-red-50 text-red-700" : "bg-green-50 text-green-700"
              }`}
            >
              {JSON.stringify(toolCall.result, null, 2)}
            </pre>
          </div>
        </div>
      )}
    </div>
  );
}
